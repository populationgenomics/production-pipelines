import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import get_batch


def cli_main():
    """
    A command-line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--dense_mt_path', help='The Dense MT to read in')
    parser.add_argument('--output_directory', help='The output directory')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(
        dense_mt_path=args.dense_mt_path,
        output_directory=args.output_directory,
    )


def main(dense_mt_path: str, output_directory: str):
    from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
    from cpg_workflows.scripts.run_pca import run

    dense_mt: hl.MatrixTable = hl.read_matrix_table(dense_mt)
    dense_mt_subset = dense_mt.head(1000)
    n_pcs = dense_mt_subset.count_cols()

    b = get_batch()
    j = dataproc_job(
        job_name='run_pca',
        function=run,
        function_path_args=dict(
            dense_mt=dense_mt_subset,
            output_directory=output_directory,
            n_pcs=n_pcs,
        ),
    )


if '__main__' == __name__:
    cli_main()
