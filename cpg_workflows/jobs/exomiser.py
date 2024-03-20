"""
jobs required for the exomiser workflow
"""
import pandas as pd
from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, authenticate_cloud_credentials_in_job, reference_path
from cpg_workflows.scripts import extract_vcf_from_mt, mt_from_vds as vds2mt, vds_from_gvcfs
from cpg_workflows.targets import Dataset, SequencingGroup
from cpg_workflows.utils import chunks, exists


def create_vds_jobs(sgids: list[SequencingGroup], out_path: str):
    """
    returns a job or None

    Args:
        sgids (list[SequencingGroup])
        out_path (str):
    """

    if exists(out_path):
        return []

    gvcf_paths = " ".join([str(s.gvcf.path) for s in sgids if s.gvcf])
    sequencing_group_names = " ".join([s.id for s in sgids])
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


def extract_mini_ped_files(
    dataset: Dataset, family_dict: dict[str, list[str]], out_prefix: Path
):
    """
    write the mini-ped for each family

    Args:
        dataset (Dataset): the dataset
        family_dict (dict[str, list[str]]): family ID to list of SG IDs
        out_prefix (Path): temp dir to write the mini-peds to
    """

    # collect all members per-family
    datasets_by_family: dict[str, list[dict]] = {}

    # limit PEDs to only those in the family_dict
    relevant_families = set(family_dict.keys())

    # query for SG entities, group by family
    for sg in dataset.get_sequencing_groups():
        ped_details = sg.pedigree.get_ped_dict()

        # skip ones we're no interested in
        if ped_details['Family.ID'] not in relevant_families:
            continue

        # add to the list for the family
        datasets_by_family.setdefault(ped_details['Family.ID'], []).append(ped_details)

    # now make the pedigree for each family
    for family, members in datasets_by_family.items():
        df = pd.DataFrame(members)
        ped_path = out_prefix / f'{family}.ped'

        # don't recreate if it exists
        if not ped_path.exists():
            with ped_path.open('w') as ped_file:
                df.to_csv(ped_file, sep='\t', index=False, header=False)


def extract_vcf_jobs(family_dict: dict[str, list[str]], mt_path: str, out_path: Path):
    """
    create the jobs to extract the VCFs from the MT

    Args:
        family_dict (dict[str, list[str]]): family ID to list of SG IDs
        mt_path (str): path to the MT
        out_path (str): temp dir to write the VCFs to

    Returns:
        list[hb.Job]: the list of jobs
    """

    vcf_jobs = []

    # find the path to the script in _this_ container
    script_path = str(extract_vcf_from_mt.__file__).removeprefix('/production-pipelines')
    for family_id, sg_ids in family_dict.items():

        vcf_target = out_path / f'{family_id}.vcf.bgz'
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
        job.command(f'python3 {script_path} {mt_path} {" ".join(sg_ids)} {vcf_target}')
        vcf_jobs.append(job)
    return vcf_jobs


def run_exomiser_batches(content_dict: dict[str, dict[str, Path]]):
    """
    run the exomiser batch
    """

    chunk_size: int = get_config()['workflow'].get('exomiser_chunk_size', 5)

    # find the exomiser references
    # first - the core ones
    clinvar_file = reference_path('exomiser_core/clinvar_whitelist')
    core_group = get_batch().read_input_group(
        **{
            'clinvar': str(clinvar_file),
            'clinvar_index': f'{clinvar_file}.tbi',
            'genome_h2': str(reference_path('exomiser_core/genome_h2')),
            'ensembl': str(reference_path('exomiser_core/ensembl_transcripts')),
            'variants': str(reference_path('exomiser_core/variants')),
        }
    )

    # cadd
    cadd_group = get_batch().read_input_group(
        **{
            'cadd_indel': str(reference_path('exomiser_cadd/indel_tsv')),
            'cadd_indel_index': str(reference_path('exomiser_cadd/indel_index')),
            'cadd_snv': str(reference_path('exomiser_cadd/snv_tsv')),
            'cadd_snv_index': str(reference_path('exomiser_cadd/snv_index')),
        }
    )

    # phenotype
    phenotype = get_batch().read_input_group(
        **{
            'pheno_db': str(reference_path('exomiser_phenotype/pheno_db')),
            'hpo_obo': str(reference_path('exomiser_phenotype/hpo_obo')),
            'rw_string': str(reference_path('exomiser_phenotype/rw_string')),
            'phenix_tar': f'{reference_path("exomiser_phenotype/phenix")}.tar.gz',
        }
    )

    # remm
    remm_group = get_batch().read_input_group(
        **{
            'remm': str(reference_path('exomiser_remm/remm_tsv')),
            'remm_index': str(reference_path('exomiser_remm/remm_tsv')),
        }
    )

    # blah
    ...
    # now chunk the jobs - load resources, then run a bunch of families
    # TODO can we load more cores and run tasks in parallel? e.g. end with &
    # UGH, we need phenopackets?!
    for instructions in chunks(content_dict, chunk_size):
        # todo ok, so we need to overwrite the application.properties file...
        # create a new job, reference the resources in a config file
        # see https://exomiser.readthedocs.io/en/latest/installation.html#linux-install
        families = list(instructions.keys())
        job = get_batch().new_bash_job(f'Run Exomiser for {families}')
        job.storage('200Gi')
        job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/exomiser:14.0.0')
        job.command(f"""
        mkdir -p data/2303_phenotype
        tar xzf {phenotype['phenix_tar']} -C data/2303_phenotype
        mv {phenotype['pheno_db']} {phenotype['hpo_obo']} {phenotype['rw_string']} data/2303_phenotype/
        ln -s {core_group} data/2303_hg38
        echo exomiser.hg38.remm-path={remm_group["remm"]} >> application.properties
        echo exomiser.hg38.cadd-snv-path={cadd_group["cadd_snv"]} >> application.properties
        echo exomiser.hg38.cadd-in-del-path={cadd_group["cadd_indel"]} >> application.properties
        tree data
        cat application.properties
        """)
