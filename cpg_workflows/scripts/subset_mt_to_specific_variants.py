#!/usr/bin/env python

"""
Take a dump TSV from seqr containing a list of variants to extract, and a MatrixTable
Creates a new subset MatrixTable to only include the variants listed in the file
"""

from argparse import ArgumentParser
from csv import DictReader

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger

logger = get_logger(__file__)


def get_variants_from_seqr(seqr_tsv: str) -> list[str]:
    """
    Read the seqr file, and return a list of variant IDs

    Args:
        seqr_tsv (str): path to the seqr file

    Returns:
        list[list[str | int]]: list of variant components
    """

    return_components = []
    with to_path(seqr_tsv).open('r') as handle:
        reader = DictReader(handle, delimiter='\t')
        for row in reader:
            return_components.append(f"chr{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}")

    return return_components


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--mt', help='Path to MatrixTable')
    parser.add_argument('--seqr', help='Path to MatrixTable')
    parser.add_argument('--output', help='Path to write MatrixTable')
    args = parser.parse_args()

    main(mt_path=args.mt, seqr_tsv=args.seqr, out_path=args.output)


def main(mt_path: str, seqr_tsv: str, out_path: str):
    """
    Subset the MatrixTable to the provided list of variants
    Args:
        mt_path (str): cohort-level matrix table from VCF.
        seqr_tsv (str): path to the seqr dump TSV
        out_path (str): where to write the results.
    """
    init_batch()

    # read the list of variants from the seqr dump TSV
    variant_list = get_variants_from_seqr(seqr_tsv)

    hl_var_list = hl.literal([hl.parse_variant(variant, reference_genome='GRCh38') for variant in variant_list])

    # read the MatrixTable
    mt = hl.read_matrix_table(mt_path)

    # subset the MatrixTable to the variants in the variant list
    mt = mt.filter_rows(hl_var_list.contains(hl.struct(locus=mt.locus, alleles=mt.alleles)))

    # write the MatrixTable
    mt.write(out_path, overwrite=True)

    logger.info(f'Subset MatrixTable written to {out_path}')


if __name__ == '__main__':
    cli_main()
