"""
input is dense mt
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch


def cli_main():

    init_batch()

    parser = ArgumentParser()
    parser.add_argument('--dense-mt', help='Path to the dense MT')
    parser.add_argument('--output-path', help='Path to the output plink files')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(dense_mt=args.dense_mt, output_path=args.output_path)


def main(dense_mt_path: str, output_path: str):

    dense_mt = hl.read_matrix_table(dense_mt_path)
    hl.export_plink(dense_mt, output_path, ind_id=dense_mt.s)


if __name__ == '__main__':
    cli_main()
