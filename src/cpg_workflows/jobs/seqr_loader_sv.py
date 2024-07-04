"""
Hail Query Batch-Backend jobs for seqr-loader.
"""

from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.query_modules import seqr_loader, seqr_loader_sv


def annotate_cohort_jobs_sv(
    vcf_path: Path,
    out_mt_path: Path,
    checkpoint_prefix: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Annotate cohort for seqr loader, SV style.
    Mostly a duplicate of the small variant version
    """

    j = get_batch().new_job('Annotate cohort', job_attrs)
    j.image(get_config()['workflow']['driver_image'])
    j.command(
        query_command(
            seqr_loader_sv,
            seqr_loader_sv.annotate_cohort_sv.__name__,
            str(vcf_path),
            str(out_mt_path),
            str(checkpoint_prefix),
            setup_gcp=True,
        ),
    )
    return j


def annotate_dataset_jobs_sv(
    mt_path: Path,
    sgids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
    exclusion_file: str | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).

    exclusion_file: str | None - this is a file path in GCP which will contain
        all SG IDs which are to be removed from consideration in this run. i.e.
        we anticipate that they will have been removed from the joint-call, and
        we will not find them, so this subsetting code will fail if it attempts
        to find them in the callset
    """
    assert sgids

    # remove excluded SG IDs
    sgids_list_path = tmp_prefix / 'sgid-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with sgids_list_path.open('w') as f:
            f.write(','.join(sgids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'

    subset_j = get_batch().new_job('subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'})
    subset_j.image(image_path('cpg_workflows'))
    subset_j.command(
        query_command(
            seqr_loader,
            seqr_loader.subset_mt_to_samples.__name__,
            str(mt_path),
            sgids,
            str(subset_mt_path),
            exclusion_file,
            setup_gcp=True,
        ),
    )

    annotate_j = get_batch().new_job('annotate dataset', (job_attrs or {}) | {'tool': 'hail query'})
    annotate_j.image(image_path('cpg_workflows'))
    annotate_j.command(
        query_command(
            seqr_loader_sv,
            seqr_loader_sv.annotate_dataset_sv.__name__,
            str(subset_mt_path),
            str(out_mt_path),
            setup_gcp=True,
        ),
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]
