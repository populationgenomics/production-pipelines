"""
Hail Query Batch-Backend jobs for seqr-loader.
"""

from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.query_modules import seqr_loader


def annotate_dataset_jobs(
    mt_path: Path,
    sequencing_group_ids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    sample_ids_list_path = tmp_prefix / 'sample-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with sample_ids_list_path.open('w') as f:
            f.write(','.join(sequencing_group_ids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'

    subset_j = get_batch().new_job('Subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'})
    subset_j.image(image_path('cpg_workflows'))
    assert sequencing_group_ids
    subset_j.command(
        query_command(
            seqr_loader,
            seqr_loader.subset_mt_to_samples.__name__,
            str(mt_path),
            sequencing_group_ids,
            str(subset_mt_path),
            setup_gcp=True,
        ),
    )
    if depends_on:
        subset_j.depends_on(*depends_on)

    annotate_j = get_batch().new_job('Annotate dataset', (job_attrs or {}) | {'tool': 'hail query'})
    annotate_j.image(image_path('cpg_workflows'))
    annotate_j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_dataset_mt.__name__,
            str(subset_mt_path),
            str(out_mt_path),
            setup_gcp=True,
        ),
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]


def cohort_to_vcf_job(
    b: Batch,
    mt_path: Path,
    out_vcf_path: Path,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None,
):
    """
    Take the single-dataset MT, and write to a VCF

    Args:
        b (hb.Batch): the batch to add jobs into
        mt_path (str): path of the AnnotateDataset MT
        out_vcf_path (Path): path to write new VCF to
        job_attrs (dict):
        depends_on (hb.Job|list[hb.Job]): jobs to depend on

    Returns:
        this single Job
    """
    from cpg_workflows.query_modules import seqr_loader

    vcf_j = b.new_job('VCF from dataset MT', (job_attrs or {}) | {'tool': 'hail query'})
    vcf_j.image(image_path('cpg_workflows'))
    vcf_j.command(
        query_command(
            seqr_loader,
            seqr_loader.vcf_from_mt_subset.__name__,
            str(mt_path),
            str(out_vcf_path),
            setup_gcp=True,
        ),
    )
    if depends_on:
        vcf_j.depends_on(*depends_on)

    return vcf_j
