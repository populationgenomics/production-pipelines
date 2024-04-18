#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

"""
Script to submit a dataproc job to run a SNP PCA

"""
import click
from cpg_utils.hail_batch import get_batch


@click.option('--vds-path', help='Path to the VDS file', required=True)
@click.option('--sample-id-file-path', help='Path to the sample id file', default=None)
@click.command()
def main(vds_path, sample_id_file_path):
    from cpg_workflows.dataproc_scripts.pca_runner import pca_runner
    from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

    dataproc_job(
        function=pca_runner,
        function_path_args=dict(vds_path=vds_path, sample_id_file_path=sample_id_file_path),
        job_name='SNP-PCA',
        num_workers=5,
    )

    get_batch(name='SNP PCA Data proc job').run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
