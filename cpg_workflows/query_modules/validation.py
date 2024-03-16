"""
Hail query method(s) for validation.
"""

import hail as hl


def single_sample_vcf_from_dataset_vcf(input_mt: str, sg_ids: list[str], out_path: str):
    """
    takes the validation datatset VCF, filters to single sample
    removes variants not relevant to this sample, and writes to VCF
    Args:
        input_mt (str): where to read the MT
        sg_ids (list[str]): one or more Sequencing Group IDs
        out_path (str): where to write the VCF to
    """
    # read full MT
    mt = hl.read_matrix_table(input_mt)

    sample_id_list_literal = hl.literal(sg_ids)

    # filter to this column
    mt = mt.filter_cols(sample_id_list_literal.contains(mt.s))

    # filter to this sample's non-ref calls
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)

    # filter out any Filter-failures
    mt = mt.filter_rows(mt.filters.length() == 0)

    hl.export_vcf(mt, out_path, tabix=True)
