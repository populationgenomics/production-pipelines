#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

"""
Script to submit a dataproc job

"""
from cpg_utils.hail_batch import get_batch
from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
from str_pca_dataproc_hail_script import pca_runner

def main():
    script = (
        f'dataproc_helper/str_pca_dataproc_hail_script.py '
        f'--file-path=gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt'
    )
    dataproc_job(
                function = pca_runner,
                function_path_args=dict(file_path='gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt'),
                job_name='STR-PCA')

    get_batch.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter