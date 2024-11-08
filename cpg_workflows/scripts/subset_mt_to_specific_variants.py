#!/usr/bin/env python

"""
Aim here is to take a list of variants, and a MatrixTable, and subset the MatrixTable
to only include the variants in the list
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger

logger = get_logger(__file__)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to MatrixTable')
    parser.add_argument('--output', help='Path to write MatrixTable')
    args = parser.parse_args()

    main(mt_path=args.input, out_path=args.output)


def main(mt_path: str, out_path: str):
    """
    Subset the MatrixTable to the provided list of variants
    Args:
        mt_path (str): cohort-level matrix table from VCF.
        out_path (str): where to write the results.
    """
    hl.init(default_reference='GRCh38')
    # init_batch()

    # read the list of variants from config, a list of strings, e.g. ['1:100000:A:T', '1:100034:A:C']
    variant_list = config_retrieve(['variant_list'])

    hl_var_list = hl.literal([
        hl.parse_variant(variant, reference_genome='GRCh38') for variant in variant_list
    ])

    # read the MatrixTable
    mt = hl.read_matrix_table(mt_path)

    # subset the MatrixTable to the variants in the variant list
    mt = mt.filter_rows(hl_var_list.contains(hl.struct(locus=mt.locus, alleles=mt.alleles)))

    # write the MatrixTable
    mt.write(out_path, overwrite=True)

    logger.info(f'Subset MatrixTable written to {out_path}')


if __name__ == '__main__':
    cli_main()