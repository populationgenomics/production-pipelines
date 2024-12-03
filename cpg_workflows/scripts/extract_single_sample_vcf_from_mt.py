"""
Extract a single sample from a MatrixTable as a VCF
Only retain filter-passing variants where this sample has a variant call
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch


def single_sample_vcf_from_dataset_vcf(
    input_mt: str,
    sample_id: str,
    out_path: str,
    clean: bool = True,
) -> None:
    """
    takes the validation datatset MatrixTable, filters to single sample
    removes variants not relevant to this sample, and writes to VCF
    Args:
        input_mt (str): where to read the MT
        sample_id (str): this Sequencing Group ID
        out_path (str): where to write the VCF to
        clean (bool): if True, remove any variant sites which are not PASS or have no calls
    """

    init_batch()

    # read full MT
    mt = hl.read_matrix_table(input_mt)

    # filter to this column
    mt = mt.filter_cols(mt.s == sample_id)

    if clean:

        # filter to this sample's non-ref calls
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

        # We need to check if this field exists on the row
        # it will for ingestion from VCF, but not for VDS derived data
        if 'filters' in mt.row:
            mt = mt.filter_rows(mt.filters.length() == 0)

        # drop the variant QC content before writing out
        mt = mt.drop('variant_qc')

    hl.export_vcf(mt, out_path, tabix=True)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='The MT to read from')
    parser.add_argument('--sample_id', help='The sample ID to filter to')
    parser.add_argument('--output', help='The VCF to write to')
    parser.add_argument(
        '--clean',
        help='If used, remove variants with no calls or not PASS',
        action='store_true',
    )
    args = parser.parse_args()

    single_sample_vcf_from_dataset_vcf(
        input_mt=args.input,
        sample_id=args.sample_id,
        out_path=args.output,
        clean=args.clean,
    )


if __name__ == '__main__':
    cli_main()
