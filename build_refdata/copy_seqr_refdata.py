"""
Copy seqr reference data to gs://cpg-reference.
"""
import os

import click
from analysis_runner import dataproc
from workflows import to_path
from workflows.batch import setup_batch

VERSION = 'v0-1'


@click.command()
def main():
    """Copy seqr reference data under ."""
    b = setup_batch(
        description='Copy seqr reference data',
        billing_project=os.environ['HAIL_BILLING_PROJECT'],
        hail_bucket=to_path('gs://cpg-reference/hail-tmp'),
    )
    dataproc.hail_dataproc_job(
        b,
        f'cpg_pipes/dataproc_scripts/copy_seqr_refdata.py {VERSION}',
        max_age='8h',
        num_secondary_workers=50,
        job_name='Copy Seqr reference data',
    )
    b.run(wait=False)


if __name__ == '__main__':
    main()
