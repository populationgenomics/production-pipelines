import logging
import os
from argparse import ArgumentParser

import hail as hl
import hailtop.batch as hb

from cpg_utils.dataproc import hail_dataproc_job

service_backend = hb.ServiceBackend(
    billing_project=os.getenv('HAIL_BILLING_PROJECT'),
    bucket=os.getenv('HAIL_BUCKET'),
)

batch = hb.Batch(name='run pca on dataproc', backend=service_backend)

hail_dataproc_job(
    batch=batch,
    script='run_pca.py',
    max_age='12h',
    num_secondary_workers=20,
    init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
    job_name='run_pca',
    worker_boot_disk_size=200,
)

batch.run(wait=False)
