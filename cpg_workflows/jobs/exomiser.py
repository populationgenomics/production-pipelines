"""
jobs required for the exomiser workflow
"""
import pandas as pd
from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, authenticate_cloud_credentials_in_job, reference_path
from cpg_workflows.scripts import extract_vcf_from_mt
from cpg_workflows.targets import Dataset
from cpg_workflows.utils import chunks


def extract_mini_ped_files(
    dataset: Dataset, family_dict: dict[str, list[str]], out_prefix: Path
):
    """
    write the mini-ped for each family

    Args:
        dataset (Dataset): the dataset
        family_dict (dict[str, list[str]]): family ID to list of SG IDs
        out_path (Path): temp dir to write the mini-peds to
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
    script_path = extract_vcf_from_mt.__file__
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
    clinvar_file = reference_path('references/exomiser_core/clinvar_whitelist')
    genome_h2 = reference_path('references/exomiser_core/genome_h2')
    ensembl = reference_path('references/exomiser_core/ensembl_transcripts')
    refseq = reference_path('references/exomiser_core/refseq_transcripts')
    ucsc = reference_path('references/exomiser_core/ucsc_transcripts')
    variants = reference_path('references/exomiser_core/variants')
    core_group = get_batch().read_input_group(
        **{
            'clinvar': str(clinvar_file),
            'clinvar_index': f'{clinvar_file}.tbi',
            'genome_h2': str(genome_h2),
            'ensembl': str(ensembl),
            'refseq': str(refseq),
            'ucsc': str(ucsc),
            'variants': str(variants),
        }
    )

    # cadd
    cadd_group = get_batch().read_input_group(
        **{
            'cadd_indel': reference_path('references/exomiser_cadd/indel_tsv'),
            'cadd_indel_index': reference_path('references/exomiser_cadd/indel_index'),
            'cadd_snv': reference_path('references/exomiser_cadd/snv_tsv'),
            'cadd_snv_index': reference_path('references/exomiser_cadd/snv_index'),
        }
    )

    # remm
    remm_group = get_batch().read_input_group(
        **{
            'remm': reference_path('references/exomiser_remm/remm'),
            'remm_index': reference_path('references/exomiser_remm/remm_tsv'),
        }
    )

    # blah
    ...
    # now chunk the jobs - load resources, then run a bunch of families
    # TODO can we load more cores and run tasks in parallel? e.g. end with &
