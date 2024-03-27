#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

"""
Script to submit a dataproc job

"""
import hail as hl
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

from str_pca_dataproc_hail_script import pca_runner


def main():

    from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

    dataproc_job(
        function=pca_runner,
        function_path_args=dict(
            file_path='gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt'
        ),
        job_name='STR-PCA',
        num_workers=5,
    )

    get_batch(name='STR PCA Data proc job').run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
