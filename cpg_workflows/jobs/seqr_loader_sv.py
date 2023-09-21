"""
Hail Query Batch-Backend jobs for seqr-loader.
"""

from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, query_command


def annotate_cohort_jobs_sv(
    b: Batch,
    vcf_path: Path,
    out_mt_path: Path,
    checkpoint_prefix: Path,
    job_attrs: dict | None = None,
    use_dataproc: bool = False,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Annotate cohort for seqr loader, SV style.
    Mostly a duplicate of the small variant version
    """
    if use_dataproc:
        from analysis_runner import dataproc

        # Script path and pyfiles should be relative to the repository root
        script = (
            'cpg_workflows/dataproc_scripts/annotate_cohort_sv.py '
            f'--vcf-path {vcf_path} '
            f'--out-mt-path {out_mt_path} '
            f'--checkpoint-prefix {checkpoint_prefix} '
        )
        pyfiles = ['cpg_workflows/query_modules']
        job_name = 'Annotate cohort'

        if cluster_id := get_config()['hail'].get('dataproc', {}).get('cluster_id'):
            # noinspection PyProtectedMember
            j = dataproc._add_submit_job(
                batch=b,
                cluster_id=cluster_id,
                script=script,
                pyfiles=pyfiles,
                job_name=job_name,
                region='australia-southeast1',
            )
        else:
            j = dataproc.hail_dataproc_job(
                b,
                script,
                max_age='24h',
                packages=[
                    'cpg_workflows',
                    'google',
                    'fsspec',
                    'gcloud',
                ],
                num_workers=2,
                num_secondary_workers=20,
                job_name=job_name,
                depends_on=depends_on,
                scopes=['cloud-platform'],
                pyfiles=pyfiles,
                init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
            )
        j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
    else:
        from cpg_workflows.query_modules import seqr_loader_sv

        j = b.new_job(f'annotate cohort', job_attrs)
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                seqr_loader_sv,
                seqr_loader_sv.annotate_cohort_sv.__name__,
                str(vcf_path),
                str(out_mt_path),
                str(checkpoint_prefix),
                setup_gcp=True,
            )
        )
        if depends_on:
            j.depends_on(*depends_on)
    return [j]


def annotate_dataset_jobs_sv(
    b: Batch,
    mt_path: Path,
    sgids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
    use_dataproc: bool = False,
    depends_on: list[Job] | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    assert sgids
    sgids_list_path = tmp_prefix / 'sgid-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with sgids_list_path.open('w') as f:
            f.write(','.join(sgids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'

    if use_dataproc:
        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        # Script path and pyfiles should be relative to the repository root
        script = (
            f'cpg_workflows/dataproc_scripts/annotate_dataset_sv.py '
            f'--mt-path {str(mt_path)} '
            f'--sgids {sgids_list_path} '
            f'--out-mt-path {out_mt_path} '
            f'--checkpoint-prefix {tmp_prefix}'
        )
        pyfiles = [
            'seqr-loading-pipelines/hail_scripts',
            'cpg_workflows/query_modules',
        ]
        job_name = 'Annotate dataset'

        if cluster_id := get_config()['hail'].get('dataproc', {}).get('cluster_id'):
            # noinspection PyProtectedMember
            j = dataproc._add_submit_job(
                batch=b,
                cluster_id=cluster_id,
                script=script,
                pyfiles=pyfiles,
                job_name=job_name,
                region='australia-southeast1',
            )
        else:
            j = dataproc.hail_dataproc_job(
                b,
                script,
                max_age='24h',
                packages=[
                    'cpg_workflows',
                    'google',
                    'fsspec',
                    'gcloud',
                ],
                num_workers=2,
                num_secondary_workers=20,
                job_name=job_name,
                depends_on=depends_on,
                scopes=['cloud-platform'],
                pyfiles=pyfiles,
                init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
            )
        j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]

    else:
        from cpg_workflows.query_modules import seqr_loader

        subset_j = b.new_job(
            f'subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'}
        )
        subset_j.image(image_path('cpg_workflows'))
        subset_j.command(
            query_command(
                seqr_loader,
                seqr_loader.subset_mt_to_samples.__name__,
                str(mt_path),
                sgids,
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
