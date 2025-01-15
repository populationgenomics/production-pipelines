"""
jobs required for the exomiser workflow
"""

import json

import pandas as pd

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path, reference_path
from cpg_utils.hail_batch import get_batch
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


def create_gvcf_to_vcf_jobs(
    proband_dict: dict[str, list[SequencingGroup]],
    previous_completions: set[str],
    out_paths: dict[str, Path],
) -> list[Job]:
    """
    Create Joint VCFs for families of SG IDs

    Args:
        proband_dict (): dict of proband ID to list of SG IDs
        previous_completions (set[str]): set of analyses we've already completed
        out_paths (): dict of family ID to output path
    Returns:
        list of Jobs
    """

    jobs: list[Job] = []

    # take each family
    for proband, members in proband_dict.items():

        # skip if already done
        if exists(out_paths[proband]) or proband in previous_completions:
            continue

        jobs.append(family_vcf_from_gvcf(members, str(out_paths[proband])))
    return jobs


def extract_mini_ped_files(proband_dict: dict[str, list[SequencingGroup]], out_paths: dict[str, Path]):
    """
    write the mini-ped for each family

    Args:
        proband_dict (dict[str, list[SequencingGroup]]): proband SG ID to list of SG IDs
        out_paths (Path): temp dir to write the mini-peds to
    """

    # query for SG entities, group by family
    for proband_id, members in proband_dict.items():

        ped_path = out_paths[proband_id]
        # don't recreate if it exists
        if not exists(ped_path):
            # make the pedigree for this family
            df = pd.DataFrame([sg.pedigree.get_ped_dict() for sg in members])
            with ped_path.open('w') as ped_file:
                df.to_csv(ped_file, sep='\t', index=False, header=False)


def make_phenopackets(proband_dict: dict[str, list[SequencingGroup]], out_paths: dict[str, Path]):
    """
    find the minimal data to run an exomiser analysis
    n.b. these are not actually phenopackets - they are a simplified version

    Args:
        proband_dict (dict[str, list[SequencingGroup]]): proband CPG ID: [family members]
        out_paths (dict[str, Path]): corresponding output paths per family
    """

    for proband, members in proband_dict.items():

        # skip if already done
        if exists(out_paths[proband]):
            continue

        # get all affected and unaffected
        affected = [sg for sg in members if str(sg.pedigree.phenotype) == '2']

        if not affected:
            raise ValueError(f'Family {proband} has no affected individuals, should not have reached here')

        # select the specific proband SG based on the ID (key in the dictionary)
        proband_sg = [sg for sg in members if sg.id == proband][0]

        hpo_term_string = proband_sg.meta['phenotypes'].get(HPO_KEY, '')

        hpo_terms = hpo_term_string.split(',')

        # https://github.com/exomiser/Exomiser/blob/master/exomiser-cli/src/test/resources/pfeiffer-family.yml
        phenopacket: dict = {'family': proband, 'proband': proband_sg.id, 'hpoIds': hpo_terms}

        with out_paths[proband].open('w') as ppk_file:
            json.dump(phenopacket, ppk_file, indent=2)


def run_exomiser(content_dict: dict[str, dict[str, Path | dict[str, Path]]]):
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
    probands = sorted(content_dict.keys())
    all_jobs = []
    for chunk_number, proband_chunk in enumerate(
        chunks(probands, config_retrieve(['workflow', 'exomiser_chunk_size'], 8)),
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

        job.command(f'echo "This job contains families {" ".join(proband_chunk)}"')

        # number of chunks should match cpu, accessible in config
        # these will all run simultaneously using backgrounded tasks and a wait
        for parallel_chunk in chunks(
            proband_chunk,
            chunk_size=config_retrieve(['workflow', 'exomiser_parallel_chunks'], 4),
        ):
            for proband in parallel_chunk:
                # read in VCF & index
                vcf = get_batch().read_input_group(
                    **{
                        f'{proband}_vcf': str(content_dict[proband]['vcf']),
                        f'{proband}_vcf_index': f'{content_dict[proband]["vcf"]}.tbi',
                    },
                )[f'{proband}_vcf']

                # read in ped & phenotype JSON
                ped = get_batch().read_input(str(content_dict[proband]['ped']))
                ppk = get_batch().read_input(str(content_dict[proband]['pheno']))

                # # this was really satisfying syntax to work out
                job.declare_resource_group(
                    **{
                        proband: {
                            'json': '{root}.json',
                            'tsv': '{root}.tsv',
                            'variants.tsv': '{root}.variants.tsv',
                            'yaml': '{root}.yaml',
                        },
                    },
                )

                # generate a config file based on the batch tmp locations
                job.command(f'python3 {exomiser_dir}/config_shuffle.py {ppk} {job[proband]["yaml"]} {ped} {vcf} ')

                # now run it, as a backgrounded process
                job.command(
                    f'java -Xmx10g -Xms4g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
                    f'--analysis {job[proband]["yaml"]} --ped {ped} '
                    f'--spring.config.location={exomiser_dir}/application.properties &',
                )

            # wait for backgrounded processes to finish, show current state
            job.command('wait && ls results')

            # move the results, then copy out
            for proband in parallel_chunk:
                job.command(f'mv results/{proband}.json {job[proband]["json"]}')
                job.command(f'mv results/{proband}.genes.tsv {job[proband]["tsv"]}')
                job.command(f'mv results/{proband}.variants.tsv {job[proband]["variants.tsv"]}')

                get_batch().write_output(
                    job[proband],
                    str(content_dict[proband]['output']).removesuffix('.tsv'),
                )
    return all_jobs
