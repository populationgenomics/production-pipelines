"""
jobs required for the exomiser workflow
"""

import json

import pandas as pd

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.scripts import collect_dataset_tsvs
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import chunks, exists

HPO_KEY: str = 'HPO Terms (present)'


def family_vcf_from_gvcf(family_members: list[SequencingGroup], out_path: str) -> Job:
    """
    Does a quick blast of gVCF -> VCF
    Strips out ref-only sites, and splits Alt, Non_Ref into just Alt
    If there are multiple members, this finishes with a merge

    Args:
        family_members (list[SequencingGroup]): member(s) to convert to VCF
        out_path (str): where to write the VCF (and implicitly, the index)

    Returns:
        the job, for dependency setting
    """

    family_ids = [sg.id for sg in family_members]
    job = get_batch().new_job(f'Generate VCF {out_path} from {family_ids}')
    job.image(image_path('bcftools'))
    if get_config()['workflow']['sequencing_type'] == 'genome':
        job.storage('10Gi')

    # read input
    family_vcfs = [get_batch().read_input(str(sg.gvcf)) for sg in family_members]

    # declare a resource group
    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # view -m 3 to strip out ref-only sites
    # norm -m -any to split Alt, Non_Ref into just Alt
    # grep -v NON_REF to remove the NON_REF sites
    # bgzip -c to write to a compressed file
    if len(family_vcfs) == 1:
        gvcf_input = family_vcfs[0]
        job.command(
            f'bcftools view -m3 {gvcf_input} | '
            f'bcftools norm -m -any | '
            f'grep -v NON_REF | '
            f'bgzip -c  > {job.output["vcf.bgz"]}',
        )
        job.command(f'tabix {job.output["vcf.bgz"]}')
        get_batch().write_output(job.output, out_path.removesuffix('.vcf.bgz'))
        return job

    # if there are multiple members, convert and merge them
    job.cpu(4)
    paths = []
    for index, gvcf_input in enumerate(family_vcfs):
        job.command(
            f'bcftools view -m3 {gvcf_input} | bcftools norm -m -any | grep -v NON_REF | bgzip -c  > {index}.vcf.bgz',
        )
        job.command(f'tabix {index}.vcf.bgz')
        paths.append(f'{index}.vcf.bgz')
    # merge the VCFs
    # -m all to merge all sites
    # -0 to replace missing with WT (potentially inaccurate, but predictable parsing in exomiser)
    # --threads 4 to use 4 threads
    # -Oz to write a compressed VCF
    job.command(f'bcftools merge {" ".join(paths)} -Oz -o {job.output["vcf.bgz"]} --threads 4 -m all -0')
    job.command(f'tabix {job.output["vcf.bgz"]}')
    get_batch().write_output(job.output, out_path.removesuffix('.vcf.bgz'))
    return job


def create_gvcf_to_vcf_jobs(families: dict[str, list[SequencingGroup]], out_paths: dict[str, Path]):
    """
    Create Joint VCFs for families of SG IDs

    Args:
        families (): dict of family ID to list of SG IDs
        out_paths (): dict of family ID to output path
    """

    jobs: list[Job] = []

    # this is the script to accomplish the same using hail VDS -> MT -> VCF
    # script_path = gvcfs_to_vcf.__file__.removeprefix('/production-pipelines')

    # take each family
    for family, members in families.items():

        # skip if already done
        if exists(out_paths[family]):
            continue

        jobs.append(family_vcf_from_gvcf(members, str(out_paths[family])))
    return jobs


def extract_mini_ped_files(family_dict: dict[str, list[SequencingGroup]], out_paths: dict[str, Path]):
    """
    write the mini-ped for each family

    Args:
        family_dict (dict[str, list[SequencingGroup]]): family ID to list of SG IDs
        out_paths (Path): temp dir to write the mini-peds to
    """

    # query for SG entities, group by family
    for family_id, members in family_dict.items():

        ped_path = out_paths[family_id]
        # don't recreate if it exists
        if not ped_path.exists():
            # make the pedigree for this family
            df = pd.DataFrame([sg.pedigree.get_ped_dict() for sg in members])
            with ped_path.open('w') as ped_file:
                df.to_csv(ped_file, sep='\t', index=False, header=False)


def make_phenopackets(family_dict: dict[str, list[SequencingGroup]], out_path: dict[str, Path]):
    """
    find the minimal data to run an exomiser analysis
    n.b. these are not actually phenopackets - they are a simplified version

    Args:
        family_dict (dict[str, list[SequencingGroup]]): members per family
        out_path (dict[str, Path]): corresponding output paths per family
    """

    for family, members in family_dict.items():

        if out_path[family].exists():
            continue

        # get all affected and unaffected
        affected = [sg for sg in members if str(sg.pedigree.phenotype) == '2']

        if not affected:
            raise ValueError(f'Family {family} has no affected individuals, should not have reached here')

        # arbitrarily select a proband for now
        proband = affected.pop()

        hpo_term_string = proband.meta['phenotypes'].get(HPO_KEY, '')

        hpo_terms = hpo_term_string.split(',')

        # https://github.com/exomiser/Exomiser/blob/master/exomiser-cli/src/test/resources/pfeiffer-family.yml
        phenopacket: dict = {'family': family, 'proband': proband.id, 'hpoIds': hpo_terms}

        with out_path[family].open('w') as ppk_file:
            json.dump(phenopacket, ppk_file, indent=2)


def run_exomiser_13(content_dict: dict[str, dict[str, Path | dict[str, Path]]]):
    """
    run the families through exomiser 13
    this is a bit of a monster, but it's a good example of how to chunk jobs
    that previous line was written by GitHub CoPilot
    This version implements Exomiser 13, with default references + CADD and REMM
    This is not recommended by the developers of Exomiser, but if we decide to
    involve CADD or REMM in the future, we'll have this code as a reference

    This was set up using a de-compressed reference bundle, which was a bit of a pain
    See the following method (Exomiser 14) for a cleaner implementation:
    - store compressed resource bundles
    - copy compressed data into image
    - inflate at runtime to build reference data folder
    """

    exomiser_version = image_path('exomiser').split(':')[-1]
    exomiser_dir = f'/exomiser/exomiser-cli-{exomiser_version}'

    # localise the exomiser references
    clinvar_file = reference_path('exomiser_core/clinvar_whitelist')
    core_group = get_batch().read_input_group(
        **{
            'clinvar': clinvar_file,
            'clinvar_index': f'{clinvar_file}.tbi',
            'genome_h2': reference_path('exomiser_core/genome_h2'),
            'ensembl': reference_path('exomiser_core/ensembl_transcripts'),
            'variants': reference_path('exomiser_core/variants'),
        },
    )

    # cadd
    cadd_group = get_batch().read_input_group(
        **{
            'cadd_indel': reference_path('exomiser_cadd/indel_tsv'),
            'cadd_indel_index': reference_path('exomiser_cadd/indel_index'),
            'cadd_snv': reference_path('exomiser_cadd/snv_tsv'),
            'cadd_snv_index': reference_path('exomiser_cadd/snv_index'),
        },
    )

    # phenotype
    phenotype = get_batch().read_input_group(
        **{
            'pheno_db': reference_path('exomiser_phenotype/pheno_db'),
            'hpo_obo': reference_path('exomiser_phenotype/hpo_obo'),
            'rw_string': reference_path('exomiser_phenotype/rw_string'),
            'phenix_tar': f'{reference_path("exomiser_phenotype/phenix")}.tar.gz',
        },
    )

    # remm
    remm_group = get_batch().read_input_group(
        **{
            'remm': reference_path('exomiser_remm/remm_tsv'),
            'remm_index': reference_path('exomiser_remm/remm_index'),
        },
    )

    # now chunk the jobs - load resources, then run a bunch of families
    families = sorted(content_dict.keys())
    all_jobs = []
    for chunk_number, family_chunk in enumerate(
        chunks(families, get_config()['workflow'].get('exomiser_chunk_size', 8)),
    ):
        # create a new job, reference the resources in a config file
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        job = get_batch().new_bash_job(f'Run Exomiser for chunk {chunk_number}')
        all_jobs.append(job)
        # higher requirement for exomiser 13 resources
        job.storage('200Gi')
        job.memory(get_config()['workflow'].get('exomiser_memory', '60Gi'))
        job.cpu(get_config()['workflow'].get('exomiser_cpu', 4))
        job.image(image_path('exomiser'))
        job.command(
            f"""
        mkdir -p {exomiser_dir}/data/2302_phenotype
        tar xzf {phenotype['phenix_tar']} -C {exomiser_dir}/data/2302_phenotype
        mv {phenotype['pheno_db']} {phenotype['hpo_obo']} {phenotype['rw_string']} {exomiser_dir}/data/2302_phenotype/
        ln -s {core_group} {exomiser_dir}/data/2302_hg38
        echo exomiser.hg38.remm-path={remm_group["remm"]} >> {exomiser_dir}/application.properties
        echo exomiser.hg38.cadd-snv-path={cadd_group["cadd_snv"]} >> {exomiser_dir}/application.properties
        echo exomiser.hg38.cadd-in-del-path={cadd_group["cadd_indel"]} >> {exomiser_dir}/application.properties
        cat {exomiser_dir}/application.properties
        """,
        )

        job.command(f'echo "This job contains families {" ".join(family_chunk)}"')

        # number of chunks should match cpu, accessible in config
        # these will all run simultaneously using backgrounded tasks and a wait
        for parallel_chunk in chunks(
            family_chunk,
            chunk_size=get_config()['workflow'].get('exomiser_parallel_chunks', 4),
        ):
            for family in parallel_chunk:
                # read in VCF & index
                vcf = get_batch().read_input_group(
                    **{
                        f'{family}_vcf': str(content_dict[family]['vcf']),
                        f'{family}_vcf_index': f'{content_dict[family]["vcf"]}.tbi',
                    },
                )[f'{family}_vcf']

                # read in ped & phenotype JSON
                ped = get_batch().read_input(str(content_dict[family]['ped']))
                ppk = get_batch().read_input(str(content_dict[family]['pheno']))

                # # this was really satisfying syntax to work out
                job.declare_resource_group(
                    **{
                        family: {
                            'json': '{root}.json',
                            'tsv': '{root}.tsv',
                            'variants.tsv': '{root}.variants.tsv',
                            'yaml': '{root}.yaml',
                        },
                    },
                )

                # generate a config file based on the batch tmp locations
                job.command(f'python3 {exomiser_dir}/config_shuffle.py {ppk} {job[family]["yaml"]} {ped} {vcf} ')

                # now run it, as a backgrounded process
                job.command(
                    f'java -Xmx10g -Xms4g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
                    f'--analysis {job[family]["yaml"]} --ped {ped} '
                    f'--spring.config.location={exomiser_dir}/application.properties &',
                )

            # wait for backgrounded processes to finish, show current state
            job.command('wait && ls results')

            # move the results, then copy out
            for family in parallel_chunk:
                job.command(f'mv results/{family}.json {job[family]["json"]}')
                job.command(f'mv results/{family}.genes.tsv {job[family]["tsv"]}')
                job.command(f'mv results/{family}.variants.tsv {job[family]["variants.tsv"]}')

                get_batch().write_output(
                    job[family],
                    str(content_dict[family]['output']).removesuffix('.tsv'),
                )
    return all_jobs


def run_exomiser_14(content_dict: dict[str, dict[str, Path | dict[str, Path]]]):
    """
    run jobs through Exomiser 14
    type hint is still wild
    """

    exomiser_version = image_path('exomiser_14').split(':')[-1]
    exomiser_dir = f'/exomiser/exomiser-cli-{exomiser_version}'

    # localise the compressed exomiser references
    inputs = get_batch().read_input_group(
        core=reference_path('exomiser_2402_core'),
        pheno=reference_path('exomiser_2402_pheno'),
    )

    # now chunk the jobs - load resources, then run a bunch of families
    families = sorted(content_dict.keys())
    all_jobs = []
    for chunk_number, family_chunk in enumerate(
        chunks(families, config_retrieve(['workflow', 'exomiser_chunk_size'], 8)),
    ):
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        job = get_batch().new_bash_job(f'Run Exomiser for chunk {chunk_number}')
        all_jobs.append(job)
        job.storage(config_retrieve(['workflow', 'exomiser_storage'], '100Gi'))
        job.memory(config_retrieve(['workflow', 'exomiser_memory'], '60Gi'))
        job.cpu(config_retrieve(['workflow', 'exomiser_cpu'], 4))
        job.image(image_path('exomiser_14'))

        # unpack references, see linux-install link above
        job.command(f'unzip {inputs}/\* -d "{exomiser_dir}/data"')

        job.command(f'echo "This job contains families {" ".join(family_chunk)}"')

        # number of chunks should match cpu, accessible in config
        # these will all run simultaneously using backgrounded tasks and a wait
        for parallel_chunk in chunks(
            family_chunk,
            chunk_size=config_retrieve(['workflow', 'exomiser_parallel_chunks'], 4),
        ):
            for family in parallel_chunk:
                # read in VCF & index
                vcf = get_batch().read_input_group(
                    **{
                        f'{family}_vcf': str(content_dict[family]['vcf']),
                        f'{family}_vcf_index': f'{content_dict[family]["vcf"]}.tbi',
                    },
                )[f'{family}_vcf']

                # read in ped & phenotype JSON
                ped = get_batch().read_input(str(content_dict[family]['ped']))
                ppk = get_batch().read_input(str(content_dict[family]['pheno']))

                # # this was really satisfying syntax to work out
                job.declare_resource_group(
                    **{
                        family: {
                            'json': '{root}.json',
                            'tsv': '{root}.tsv',
                            'variants.tsv': '{root}.variants.tsv',
                            'yaml': '{root}.yaml',
                        },
                    },
                )

                # generate a config file based on the batch tmp locations
                job.command(f'python3 {exomiser_dir}/config_shuffle.py {ppk} {job[family]["yaml"]} {ped} {vcf} ')

                # now run it, as a backgrounded process
                job.command(
                    f'java -Xmx10g -Xms4g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
                    f'--analysis {job[family]["yaml"]} --ped {ped} '
                    f'--spring.config.location={exomiser_dir}/application.properties &',
                )

            # wait for backgrounded processes to finish, show current state
            job.command('wait && ls results')

            # move the results, then copy out
            for family in parallel_chunk:
                job.command(f'mv results/{family}.json {job[family]["json"]}')
                job.command(f'mv results/{family}.genes.tsv {job[family]["tsv"]}')
                job.command(f'mv results/{family}.variants.tsv {job[family]["variants.tsv"]}')

                get_batch().write_output(
                    job[family],
                    str(content_dict[family]['output']).removesuffix('.tsv'),
                )
    return all_jobs


def generate_seqr_summary(tsv_dict: dict[str, Path], project: str, output: str):
    """
    Generate the summary TSV for Seqr by combining all per-family TSVs

    Args:
        tsv_dict (dict[str, Path]):
        project ():
        output ():
    """

    # path to the python script
    script_path = collect_dataset_tsvs.__file__.removeprefix('/production-pipelines')

    # read all the input files in
    files_read_in = [get_batch().read_input(str(tsv)) for tsv in tsv_dict.values()]

    job = get_batch().new_bash_job(f'Aggregate TSVs for {project}')
    job.storage('10Gi')
    job.image(get_config()['workflow']['driver_image'])
    job.command(f'python3 {script_path} --project {project} --input {" ".join(files_read_in)} --output {job.output} ')
    get_batch().write_output(job.output, output)
    return [job]
