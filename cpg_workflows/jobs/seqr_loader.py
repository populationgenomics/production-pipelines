"""
Hail Query Batch-Backend jobs for seqr-loader.
"""

from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, query_command


def annotate_cohort_jobs(
    b: Batch,
    vcf_path: Path,
    out_mt_path: Path,
    checkpoint_prefix: Path,
    vep_ht_path: Path,
    siteonly_vqsr_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
    use_dataproc: bool = True,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Annotate cohort for seqr loader.
    """
    if use_dataproc:
        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        script_path = (
            to_path(__file__).parent / 'dataproc_scripts' / 'annotate_cohort.py'
        )

        j = dataproc.hail_dataproc_job(
            b,
            f'{script_path} '
            f'--vcf-path {vcf_path} '
            f'--vep-ht-path {vep_ht_path} '
            f'--out-mt-path {out_mt_path} '
            f'--checkpoint-prefix {checkpoint_prefix} '
            + (
                f'--siteonly-vqsr-vcf-path {siteonly_vqsr_vcf_path} '
                if siteonly_vqsr_vcf_path
                else ''
            ),
            max_age='24h',
            packages=[
                'cpg_workflows',
                'google',
                'fsspec',
                'gcloud',
            ],
            num_workers=2,
            num_secondary_workers=20,
            job_name=f'Annotate cohort',
            depends_on=depends_on,
            scopes=['cloud-platform'],
            pyfiles=[
                'seqr-loading-pipelines/hail_scripts',
                'cpg_workflows/query_modules',
            ],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        )
        j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
    else:
        from cpg_workflows.query_modules import seqr_loader

        j = b.new_job(f'annotate cohort', job_attrs)
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                seqr_loader,
                seqr_loader.annotate_cohort.__name__,
                str(vcf_path),
                str(out_mt_path),
                str(vep_ht_path),
                str(siteonly_vqsr_vcf_path),
                str(checkpoint_prefix),
                setup_gcp=True,
            )
        )
        if depends_on:
            j.depends_on(*depends_on)
    return [j]


def annotate_dataset_jobs(
    b: Batch,
    mt_path: Path,
    sample_ids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
    use_dataproc: bool = True,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    sample_ids_list_path = tmp_prefix / 'sample-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with sample_ids_list_path.open('w') as f:
            f.write(','.join(sample_ids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'

    if use_dataproc:
        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        script_path = (
            to_path(__file__).parent / 'dataproc_scripts' / 'annotate_dataset.py'
        )

        j = dataproc.hail_dataproc_job(
            b,
            f'{script_path} '
            f'--mt-path {mt_path} '
            f'--sample-ids {sample_ids_list_path} '
            f'--out-mt-path {out_mt_path} '
            f'--checkpoint-prefix {tmp_prefix}',
            max_age='24h',
            packages=[
                'cpg_workflows',
                'google',
                'fsspec',
                'gcloud',
            ],
            num_workers=2,
            num_secondary_workers=20,
            job_name=f'Annotate dataset',
            depends_on=depends_on,
            scopes=['cloud-platform'],
            pyfiles=[
                'seqr-loading-pipelines/hail_scripts',
                'query_modules',
            ],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        )
        j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]

    else:
        from cpg_workflows.query_modules import seqr_loader

        subset_j = b.new_job(
            f'subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'}
        )
        subset_j.image(image_path('cpg_workflows'))
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
        if depends_on:
            subset_j.depends_on(*depends_on)

        annotate_j = b.new_job(
            f'annotate dataset', (job_attrs or {}) | {'tool': 'hail query'}
        )
        annotate_j.image(image_path('cpg_workflows'))
        annotate_j.command(
            query_command(
                seqr_loader,
                seqr_loader.annotate_dataset_mt.__name__,
                str(subset_mt_path),
                str(out_mt_path),
                str(tmp_prefix),
                setup_gcp=True,
            )
        )
        annotate_j.depends_on(subset_j)
        jobs = [subset_j, annotate_j]

    return jobs
