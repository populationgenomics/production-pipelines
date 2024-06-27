#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

"""
Script to submit a dataproc job

"""
from cpg_utils.hail_batch import get_batch


def main():
    from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
    from cpg_workflows.large_cohort.str_pca_dataproc_hail_script import pca_runner

    dataproc_job(
        function=pca_runner,
        function_path_args=dict(
            file_path='gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt',
        ),
        job_name='STR-PCA',
        num_workers=20,
    )

    get_batch(name='STR PCA Data proc job').run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
