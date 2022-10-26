"""
Hail Query Batch-Backend jobs for seqr-loader.
"""

from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import Path
from cpg_utils.hail_batch import image_path, genome_build, query_command

from cpg_workflows.query_modules import seqr_loader


def annotate_cohort_jobs(
    b: Batch,
    vcf_path: Path,
    siteonly_vqsr_vcf_path: Path,
    vep_ht_path: Path,
    out_mt_path: Path,
    checkpoint_prefix: Path,
    sequencing_type: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Annotate cohort for seqr loader.
    """
    j = b.new_job(f'annotate cohort', job_attrs)
    j.image(image_path('cpg_utils'))
    j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(vcf_path),
            str(siteonly_vqsr_vcf_path),
            str(vep_ht_path),
            str(out_mt_path),
            overwrite,
            genome_build(),
            sequencing_type,
            str(checkpoint_prefix),
            setup_gcp=True,
        )
    )
    return [j]


def annotate_dataset_jobs(
    b: Batch,
    mt_path: Path,
    sample_ids: list[str],
    output_mt_path: Path,
    tmp_bucket: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    subset_mt_path = tmp_bucket / 'cohort-subset.mt'
    subset_j = b.new_job(
        f'subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'}
    )
    subset_j.image(image_path('cpg_utils'))
    assert sample_ids
    subset_j.command(
        query_command(
            seqr_loader,
            seqr_loader.subset_mt_to_samples.__name__,
            str(mt_path),
            sample_ids,
            str(subset_mt_path),
            setup_gcp=True,
        )
    )

    annotate_j = b.new_job(
        f'annotate dataset', (job_attrs or {}) | {'tool': 'hail query'}
    )
    annotate_j.image(image_path('cpg_utils'))
    annotate_j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_dataset_mt.__name__,
            str(subset_mt_path),
            str(output_mt_path),
            str(tmp_bucket),
            overwrite,
            setup_gcp=True,
        )
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]
