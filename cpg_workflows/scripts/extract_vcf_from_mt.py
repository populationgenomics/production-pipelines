#! /usr/bin/env python3

"""
Takes a MatrixTable, a list of sample IDs, and an output path
Writes a VCF representation containing only the requested samples
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger


def extract_vcf_from_dataset_vcf(input_mt: str, sg_ids: list[str], out_path: str):
    """
    takes a datatset VCF, filters to sample subset
    removes variants not relevant to subset, and writes to VCF
    the VCF write is direct to GCP

    Args:
        input_mt (str): input MT
        sg_ids (list[str]): one or more Sequencing Group IDs
        out_path (str): where to write the VCF to
    """

    # read full MT
    mt = hl.read_matrix_table(input_mt)

    samples_in_mt = set(mt.s.collect())
    requested_samples = set(sg_ids)
    if missing_samples := (requested_samples - samples_in_mt):
        raise ValueError(f'Requested samples {missing_samples} not found in MT')

    # change the string list to a hail object version
    sample_id_list_literal = hl.literal(sg_ids)

    # filter to this column
    mt = mt.filter_cols(sample_id_list_literal.contains(mt.s))

    # filter to this sample's non-ref calls
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    # drop gvcf_info - a dict, can't sit in FORMAT
    # only relevant when doing gVCF -> VDS -> MT
    if 'gvcf_info' in mt.row:
        mt = mt.drop('gvcf_info')
    mt = mt.drop('variant_qc')

    hl.export_vcf(mt, out_path, tabix=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('input_mt', help='Input MatrixTable')
    parser.add_argument('sg_ids', nargs='+', help='One or more Sequencing Group IDs')
    parser.add_argument('out_path', help='Output VCF path')
    args = parser.parse_args()

    get_logger(__file__).info(f'Extracting VCF from {args.input_mt} for {args.sg_ids}, writing to {args.out_path}')

    init_batch()
    extract_vcf_from_dataset_vcf(args.input_mt, args.sg_ids, args.out_path)
