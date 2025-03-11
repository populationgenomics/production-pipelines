#! /usr/bin/env python3

"""
Takes a MatrixTable and two output paths
Writes two representations - full VCF and a sites-only version

This script takes a BED file and a MatrixTable as input,
The output contains only the variants that overlap with the BED file
All existing info fields are dropped, and replaced with the callset
AC / AN / AF

This removes any Filtered variants, due to an issue between the header and content
 - when we apply VQSR annotations we pull in Filters, but we don't pull the corresponding header lines
 - this means that the MT contains rows which aren't explained in the header, causing some tools to fail
 - for now it's easier to just remove these rows - reconsider if we use this properly
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import genome_build
from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger


def extract_vcf_from_mt(
    mt_path: str,
    out: str,
    sites_only: str,
    bed: str,
) -> None:
    """
    takes a datatset MT
    Limits the MT to the regions specified in the BED file

    Args:
        mt_path (str): input MT
        out (str): Full VCF destination
        sites_only (str): Sites-only VCF destination
        bed (str): BED file containing regions to keep
    """

    init_batch()

    # remote-read of the BED file, skipping any contigs not in the reference genome
    # the Ensembl data wasn't filtered to remove non-standard contigs
    limited_region = hl.import_bed(bed, reference_genome=genome_build(), skip_invalid_intervals=True)

    # read full MT
    mt = hl.read_matrix_table(mt_path)

    # filter to PASS only (this is a temporary measure until we insert VQSR headers)
    mt = mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))

    # filter to overlaps with the BED file
    mt = mt.filter_rows(hl.is_defined(limited_region[mt.locus]))

    # replace the existing INFO block to just have AC/AN/AF - no other carry-over
    # this is based on the structure we already achieved in annotate_cohort
    mt = mt.annotate_rows(
        info=hl.struct(
            AF=mt.AF,
            AN=mt.AN,
            AC=mt.AC,
        ),
    )

    mt.describe()

    # write a full VCF
    hl.export_vcf(mt, out, tabix=True)
    hl.export_vcf(mt.rows(), sites_only, tabix=True)


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--mt', help='Input MatrixTable', required=True)
    parser.add_argument('--out', help='Path to full output VCF', required=True)
    parser.add_argument('--sites_only', help='Path to output sites-only VCF', required=True)
    parser.add_argument('--bed', help='Region BED file', required=True)
    args = parser.parse_args()

    get_logger(__file__).info(f'Extracting VCF from {args.mt}, limited to {args.bed}')
    get_logger(__file__).info(f'Writing full VCF to {args.out}')
    get_logger(__file__).info(f'Writing sites-only VCF to {args.sites_only}')

    extract_vcf_from_mt(
        mt_path=args.mt,
        out=args.out,
        sites_only=args.sites_only,
        bed=args.bed,
    )


if __name__ == '__main__':
    cli_main()
