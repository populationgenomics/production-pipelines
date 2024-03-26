"""
jobs required for the exomiser workflow
"""

import json

import pandas as pd

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, image_path, reference_path
from cpg_workflows.scripts import extract_vcf_from_mt, vds_from_gvcfs
from cpg_workflows.scripts import mt_from_vds as vds2mt
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import chunks, exists

HPO_KEY: str = 'HPO Terms (present)'


def create_vds_jobs(sgids: list[SequencingGroup], out_path: str):
    """
    Create VDS from a number of SG IDs
    At larger cohort sizes this won't work, as the number of
    arguments may hit a threshold

    Mitigate by writing a temp file and reading into script container

    Args:
        sgids (list[SequencingGroup])
        out_path (str):
    """

    if exists(out_path):
        return []

    gvcf_paths = " ".join([str(s.gvcf.path) for s in sgids if s.gvcf])
    sequencing_group_names = " ".join([s.id for s in sgids])

    # need a sexier way of doing this... currently there's a disconnect
    # between how the cpg_workflows image is built, compared with
    # how this driver image git checkout is done, meaning a different filepath
    vds_script = str(vds_from_gvcfs.__file__).removeprefix('/production-pipelines')

    job = get_batch().new_job('Make VDS from gVCFs')
    authenticate_cloud_credentials_in_job(job)
    job.image(get_config()['workflow']['driver_image'])
    job.command(f'python3 {vds_script} --gvcfs {gvcf_paths} --sgids {sequencing_group_names} --out {out_path}')
    job.storage('10Gi')
    job.cpu(4)
    job.memory('standard')
    return job


def mt_from_vds(vds_path: str, out_path: str):
    """
    create the MT from the VDS

    Args:
        vds_path (str): path to the VDS
        out_path (str): path to the MT
    """

    if exists(out_path):
        return []

    # find the path to the script in _this_ container
    script_path = str(vds2mt.__file__).removeprefix('/production-pipelines')
    job = get_batch().new_job('Make MT from VDS')
    job.image(get_config()['workflow']['driver_image'])
    job.command(f'python3 {script_path} --vds {vds_path} --mt {out_path}')
    job.storage('10Gi')
    job.cpu(4)
    job.memory('standard')
    return job


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


def extract_vcf_jobs(family_dict: dict[str, list[SequencingGroup]], mt_path: str, out_path: dict[str, Path]):
    """
    create the jobs to extract the VCFs from the MT

    Args:
        family_dict (dict[str, list[str]]): family ID to list of SG IDs
        mt_path (str): path to the MT
        out_path (dict): path to each family's output file

    Returns:
        list[hb.Job]: the list of jobs
    """

    vcf_jobs = []

    # find the path to the script in _this_ container
    script_path = str(extract_vcf_from_mt.__file__).removeprefix('/production-pipelines')
    for family_id, sg_ids in family_dict.items():

        vcf_target = out_path[family_id]
        if vcf_target.exists():
            continue

        job = get_batch().new_bash_job(f'Extract VCF for {family_id}')

        # set _this_ image as the execution environment
        job.image(get_config()['workflow']['driver_image'])

        # no need for this - hail will write directly
        # job.declare_resource_group(
        #     output={
        #         'vcf.bgz': '{root}.vcf.bgz',
        #         'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        #     }
        # )

        authenticate_cloud_credentials_in_job(job)
        sgid_ids = ' '.join([sg.id for sg in family_dict[family_id]])
        job.command(f'python3 {script_path} {mt_path} {sgid_ids} {vcf_target}')
        vcf_jobs.append(job)
    return vcf_jobs


def make_phenopackets(family_dict: dict[str, list[SequencingGroup]], out_path: dict[str, Path]):
    """
    find the minimal data to run an exomiser analysis
    n.b. these are not actually phenopackets - they are a simplified version

    Args:
        family_dict (dict[str, list[SequencingGroup]]): members per family
        out_path (dict[str, Path]): corresponding output paths per family
    """

    for family, members in family_dict.items():

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


def run_exomiser_batches(content_dict: dict[str, dict[str, Path]]):
    """
    run the exomiser batch
    """

    exomiser_version = image_path('exomiser').split(':')[-1]
    exomiser_dir = f'/exomiser/exomiser-cli-{exomiser_version}'

    # localise the exomiser references
    clinvar_file = reference_path('exomiser_core/clinvar_whitelist')
    core_group = get_batch().read_input_group(
        **{
            'clinvar': str(clinvar_file),
            'clinvar_index': f'{clinvar_file}.tbi',
            'genome_h2': str(reference_path('exomiser_core/genome_h2')),
            'ensembl': str(reference_path('exomiser_core/ensembl_transcripts')),
            'variants': str(reference_path('exomiser_core/variants')),
        },
    )

    # cadd
    cadd_group = get_batch().read_input_group(
        **{
            'cadd_indel': str(reference_path('exomiser_cadd/indel_tsv')),
            'cadd_indel_index': str(reference_path('exomiser_cadd/indel_index')),
            'cadd_snv': str(reference_path('exomiser_cadd/snv_tsv')),
            'cadd_snv_index': str(reference_path('exomiser_cadd/snv_index')),
        },
    )

    # phenotype
    phenotype = get_batch().read_input_group(
        **{
            'pheno_db': str(reference_path('exomiser_phenotype/pheno_db')),
            'hpo_obo': str(reference_path('exomiser_phenotype/hpo_obo')),
            'rw_string': str(reference_path('exomiser_phenotype/rw_string')),
            'phenix_tar': f'{reference_path("exomiser_phenotype/phenix")}.tar.gz',
        },
    )

    # remm
    remm_group = get_batch().read_input_group(
        **{
            'remm': str(reference_path('exomiser_remm/remm_tsv')),
            'remm_index': str(reference_path('exomiser_remm/remm_index')),
        },
    )

    # now chunk the jobs - load resources, then run a bunch of families
    families = sorted(content_dict.keys())
    all_jobs = []
    for family_chunk in chunks(families, get_config()['workflow'].get('exomiser_chunk_size', 8)):
        # create a new job, reference the resources in a config file
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        job = get_batch().new_bash_job(f'Run Exomiser for {family_chunk}')
        all_jobs.append(job)
        job.storage(get_config()['workflow'].get('exomiser_storage', '200Gi'))
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
        set -x
        tree -l data
        """,
        )

        # run the example data
        # job.command(
        #     f'java -Xmx10g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
        #     f'--analysis examples/test-analysis-multisample.yml && '
        #     'ls results && '
        #     'cat results/Pfeiffer_exomiser.json && '
        #     'cat results/Pfeiffer-quartet-hiphive-exome-PASS_ONLY.json'
        # )

        for parallel_chunk in chunks(family_chunk, chunk_size=get_config()['workflow'].get('exomiser_chunk_size', 8)):
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
                job.declare_resource_group(**{family: {'json': '{root}.json', 'yaml': '{root}.yaml'}})

                # generate a config file based on the batch tmp locations
                job.command(f'python3 {exomiser_dir}/config_shuffle.py {ppk} {job[family]["yaml"]} {ped} {vcf} ')
                job.command(f'cat {job[family]["yaml"]}')

                # now run it, as a backgrounded process
                job.command(
                    f'java -Xmx10g -Xms4g -jar {exomiser_dir}/exomiser-cli-{exomiser_version}.jar '
                    f'--analysis {job[family]["yaml"]} --ped {ped} '
                    f'--spring.config.location={exomiser_dir}/application.properties &',
                )

            # wait for backgrounded processes to finish, show current state
            job.command('wait && ls *')

            # move the results, then copy out
            # the output-prefix value can't take a path with a / in it, so we can't use the resource group
            for family in parallel_chunk:
                job.command(f'mv results/{family}.json {job[family]["json"]}')
                get_batch().write_output(
                    job[family],
                    str(content_dict[family]['output']).removesuffix('.json'),
                )
    return all_jobs
