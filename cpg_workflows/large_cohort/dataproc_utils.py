"""
Utils for submitting the workflow to a Dataproc cluster.
"""

import math
from analysis_runner import dataproc
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config

from cpg_workflows.batch import get_batch


DATAPROC_PACKAGES = [
    'cpg-utils',
    'click',
    'gnomad',
    'google',
    'slackclient',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]


# GCP quota on number of cores
MAX_PRIMARY_WORKERS = 50


def dataproc_job(
    job_name: str,
    function,
    function_path_args: dict[str, Path],
    function_str_args: list[str] | None = None,
    preemptible: bool = True,
    num_workers: int | None = None,
    depends_on: list[Job | None] = None,
    autoscaling_policy: str | None = None,
    long: bool = False,
    worker_boot_disk_size: int | None = None,
    secondary_worker_boot_disk_size: int | None = None,
) -> Job:
    """
    Submit script as a dataproc job.
    """
    from cpg_workflows import large_cohort
    from cpg_workflows.large_cohort import dataproc_script

    package_path = to_path(large_cohort.__file__).parent
    script_path = to_path(dataproc_script.__file__)
    rel_script_path = script_path.relative_to(package_path.parent)

    script = (
        f'{rel_script_path} '
        f'{function.__module__} {function.__name__} '
        f'{" ".join([f"-p {p}" for p in function_path_args.values()])} '
        f'{" ".join([a for a in function_str_args or []])} '
    )

    if cluster_id := get_config()['hail'].get('dataproc', {}).get('cluster_id'):
        # noinspection PyProtectedMember
        return dataproc._add_submit_job(
            batch=get_batch(),
            cluster_id=cluster_id,
            script=script,
            pyfiles=[large_cohort.__name__],
            job_name=job_name,
            region='australia-southeast1',
        )

    if num_workers is None:
        num_workers = get_config()['workflow']['scatter_count']

    max_age = '24h' if not long else '48h'

    depends_on = depends_on or []
    depends_on = [j for j in depends_on if j is not None]

    if autoscaling_policy:
        num_secondary_workers = 0
        num_primary_workers = 0
    elif preemptible:
        num_secondary_workers = num_workers
        # number of primary workers has to be 5-10% of the number of secondary workers:
        # see Tim's comment in https://discuss.hail.is/t/all-nodes-are-unhealthy/1764
        # using 8% here:
        num_primary_workers = int(math.ceil(num_secondary_workers * 0.08))
    else:
        num_secondary_workers = 0
        num_primary_workers = num_workers

    use_highmem_workers = get_config()['workflow'].get('highmem_workers')

    # 2 is the minimal number of primary workers for dataproc cluster:
    num_primary_workers = max(num_primary_workers, 2)

    # limiting the number of primary workers to avoid hitting GCP quota:
    num_primary_workers = min(num_primary_workers, MAX_PRIMARY_WORKERS)

    return dataproc.hail_dataproc_job(
        get_batch(),
        script=script,
        job_name=job_name,
        max_age=max_age,
        packages=DATAPROC_PACKAGES,
        num_secondary_workers=num_secondary_workers,
        num_workers=num_primary_workers,
        autoscaling_policy=autoscaling_policy,
        depends_on=depends_on,
        worker_machine_type='n1-highmem-8' if use_highmem_workers else 'n1-standard-8',
        worker_boot_disk_size=worker_boot_disk_size,
        secondary_worker_boot_disk_size=secondary_worker_boot_disk_size,
        pyfiles=['large_cohort'],
        init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    )
