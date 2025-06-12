#!/usr/bin/env python

"""
Convert to site-only table and annotate with AS fields.
"""
import logging

import hail as hl

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, dataset_path, output_path
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import can_reuse
from gnomad.utils.sparse_mt import default_compute_info
from gnomad.utils.vcf import adjust_vcf_incompatible_types


def extract_as_pls(
    lpl_expr: hl.expr.ArrayExpression,
    allele_idx: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Extract PLs for a specific allele from an LPL array expression.

    PL/LPL represents the normalized Phred-scaled likelihoods of the possible
    genotypes from all considered alleles (or local alleles).

    If three alleles are considered, LPL genotype indexes are:
    [0/0, 0/1, 1/1, 0/2, 1/2, 2/2].

    If we want to extract the PLs for each alternate allele, we need to extract:
        - allele 1: [0/0, 0/1, 1/1]
        - allele 2: [0/0, 0/2, 2/2]

    Example:
        - LPL: [138, 98, 154, 26, 0, 14]
        - Extract allele 1 PLs: [0/0, 0/1, 1/1] -> [138, 98, 154]
        - Extract allele 2 PLs: [0/0, 0/2, 2/2] -> [138, 26, 14]

    :param lpl_expr: LPL ArrayExpression.
    :param allele_idx: The index of the alternate allele to extract PLs for.
    :return: ArrayExpression of PLs for the specified allele.
    """
    calls_to_keep = hl.array(
        [hl.call(0, 0), hl.call(0, allele_idx), hl.call(allele_idx, allele_idx)],
    )
    return calls_to_keep.map(lambda c: lpl_expr[c.unphased_diploid_gt_index()])


def recompute_as_qualapprox_from_lpl(mt: hl.MatrixTable) -> hl.expr.ArrayExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Recompute AS_QUALapprox from LPL.

    QUALapprox is the (Phred-scaled) probability that all reads at the site are hom-ref,
    so QUALapprox is PL[0]. To get the QUALapprox for just one allele, pull out the
    PLs for just that allele, then normalize by subtracting the smallest element from
    all the entries (so the best genotype is 0) and then use the normalized PL[0]
    value for that allele's QUALapprox.

    .. note::

        - The first element of AS_QUALapprox is always None.
        - If the allele is a star allele, we set QUALapprox for that allele to 0.
        - If GQ == 0 and PL[0] for the allele == 1, we set QUALapprox for the allele
          to 0.

    Example:
        Starting Values:
            - alleles: [‘G’, ‘*’, ‘A’, ‘C’, ‘GCTT’, ‘GT’, ‘T’]
            - LGT: 1/2
            - LA: [0, 1, 6]
            - LPL: [138, 98, 154, 26, 0, 14]
            - QUALapprox: 138

        Use `extract_as_pls` to get PLs for each allele:
            - allele 1: [138, 98, 154]
            - allele 2: [138, 26, 14]

        Normalize PLs by subtracting the smallest element from all the PLs:
            - allele 1: [138-98, 98-98, 154-98] -> [40, 0, 56]
            - allele 2: [138-14, 26-14, 14-14] -> [124, 12, 0]

        Use the first element of the allele specific PLs to generate AS_QUALapprox:
        [None, 40, 124]

        Set QUALapprox to 0 for the star allele: [None, 0, 124]

    :param mt: Input MatrixTable.
    :return: AS_QUALapprox ArrayExpression recomputed from LPL.
    """
    return hl.enumerate(mt.LA).map(
        lambda i: (
            hl.case()
            .when(mt.alleles[i[1]] == "*", 0)
            .when(
                i[0] > 0,
                hl.bind(
                    lambda pl_0: hl.if_else((mt.GQ == 0) & (pl_0 == 1), 0, pl_0),
                    hl.bind(lambda x: x[0] - hl.min(x), extract_as_pls(mt.LPL, i[0])),
                ),
            )
            .or_missing()
        ),
    )


def correct_as_annotations(
    mt: hl.MatrixTable,
    set_to_missing: bool = False,
) -> hl.expr.StructExpression:
    """
    Function pulled from `gnomad_qc/v4/annotations/generate_variant_qc_annotations.py`.

    Correct allele specific annotations that are longer than the length of LA.

    For some entries in the MatrixTable, the following annotations are longer than LA,
    when they should be the same length as LA:

        - AS_SB_TABLE
        - AS_RAW_MQ
        - AS_RAW_ReadPosRankSum
        - AS_RAW_MQRankSum

    This function corrects these annotations by either dropping the alternate allele
    with the index corresponding to the min value of AS_RAW_MQ, or setting them to
    missing if `set_to_missing` is True.

    :param mt: Input MatrixTable.
    :param set_to_missing: Whether to set the annotations to missing instead of
        correcting them.
    :return: StructExpression with corrected allele specific annotations.
    """
    annotations_to_correct = [
        "AS_SB_TABLE",
        "AS_RAW_MQ",
        "AS_RAW_ReadPosRankSum",
        "AS_RAW_MQRankSum",
    ]
    annotations_to_correct_dict = {a: mt.gvcf_info[a] for a in annotations_to_correct}

    # Identify index corresponding to min of AS_RAW_MQ, skipping the reference allele.
    as_raw_mq_no_ref = mt.gvcf_info.AS_RAW_MQ[1:]
    idx_remove = as_raw_mq_no_ref.index(hl.min(as_raw_mq_no_ref)) + 1

    corrected_annotations = {
        a: hl.if_else(
            hl.len(expr) > hl.len(mt.LA),
            hl.or_missing(
                not set_to_missing,
                expr[:idx_remove].extend(expr[idx_remove + 1 :]),
            ),
            expr,
        )
        for a, expr in annotations_to_correct_dict.items()
    }

    return hl.struct(**corrected_annotations)


def run(
    vds_path: str,
    sample_qc_ht_path: str,
    relateds_to_drop_ht_path: str,
    out_as_vcf_path: str,
    out_quasi_vcf_path: str,
    out_ht_path: str,
    out_ht_pre_vcf_adjusted_path: str,
):
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    vds = hl.vds.read_vds(str(vds_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))

    as_info_ht, quasi_info_ht = vds_to_site_only_ht(
        vds=vds,
        sample_qc_ht=sample_qc_ht,
        relateds_to_drop_ht=relateds_to_drop_ht,
        out_ht_path=out_ht_path,
        out_ht_pre_vcf_adjusted_path=out_ht_pre_vcf_adjusted_path,
    )
    logging.info(f'Writing AS VCF to {out_as_vcf_path} and quasi-AS VCF to {out_quasi_vcf_path}')
    assert to_path(out_as_vcf_path).suffix == '.bgz'
    assert to_path(out_quasi_vcf_path).suffix == '.bgz'
    hl.export_vcf(as_info_ht, str(out_as_vcf_path), tabix=True)
    hl.export_vcf(quasi_info_ht, str(out_quasi_vcf_path), tabix=True)


def build_info_ht(ht: hl.Table, extra_field: str) -> hl.Table:
    ht = ht.select('AC_info', extra_field, 'site_info', 'lowqual', 'AS_lowqual')
    ht = ht.annotate(
        info=hl.struct(**ht.AC_info, **ht[extra_field], **ht.site_info),
    )
    ht = ht.drop(
        'AC_info',
        extra_field,
        'site_info',
    )

    ht = adjust_vcf_incompatible_types(
        ht,
        # With default INFO_VCF_AS_PIPE_DELIMITED_FIELDS, AS_VarDP will be converted
        # into a pipe-delimited value e.g.: VarDP=|132.1|140.2
        # which breaks VQSR parser (it doesn't recognise the delimiter and treats
        # it as an array with a single string value "|132.1|140.2", leading to
        # an IndexOutOfBound exception when trying to access value for second allele)
        pipe_delimited_annotations=[],
    )
    return ht


def vds_to_site_only_ht(
    vds: hl.vds.VariantDataset,
    sample_qc_ht: hl.Table,
    relateds_to_drop_ht: hl.Table,
    out_ht_path: str,
    out_ht_pre_vcf_adjusted_path: str,
) -> hl.Table:
    """
    Convert VDS into sites-only VCF-ready table.
    """
    if can_reuse(out_ht_path):
        return hl.read_table(str(out_ht_path))

    # Passing un-densified MatrixTable to default_compute_info
    mt = vds.variant_data

    mt = mt.filter_cols(hl.len(sample_qc_ht[mt.col_key].filters) > 0, keep=False)
    mt = mt.filter_cols(hl.is_defined(relateds_to_drop_ht[mt.col_key]), keep=False)

    # Compute and checkpoint the allele specific info annotations after recomputing
    # AS_QUALapprox from LPL, and fixing the length of AS_SB_TABLE, AS_RAW_MQ,
    # AS_RAW_ReadPosRankSum and AS_RAW_MQRankSum.
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))
    correct_mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.annotate(
            AS_QUALapprox=recompute_as_qualapprox_from_lpl(mt),
            **correct_as_annotations(mt),
        ),
    )

    correct_mt = _filter_rows_and_add_tags(correct_mt)
    # TODO: Figure out the right output pathing function for this checkpoint
    correct_mt = correct_mt.checkpoint(
        dataset_path('MakeSiteOnly/corrected.mt', category='tmp'),
        overwrite=True,
    )
    info_ht: hl.Table = _create_info_ht(correct_mt, n_partitions=mt.n_partitions())
    info_ht = info_ht.checkpoint(out_ht_pre_vcf_adjusted_path, overwrite=True)

    # AS_info and site_info are stored in the info_ht.info struct.
    # Separate info_ht.info fields into AC_info, AS_info, quasi_info, and site_info.
    AS_info = info_ht.info
    AS_keys = hl.eval(
        hl.array(list(AS_info.keys())).group_by(
            lambda x: (hl.switch(x[:2]).when("AS", "AS_info").when("AC", "AC_info").default("site_info")),
        ),
    )
    AS_info = {k: AS_info.select(*v) for k, v in AS_keys.items()}

    info_ht = info_ht.select(
        AC_info=AS_info["AC_info"],
        AS_info=AS_info["AS_info"],
        site_info=AS_info["site_info"],
        quasi_info=info_ht["quasi_info"],
        lowqual=info_ht.lowqual,
        AS_lowqual=info_ht.AS_lowqual,
    )

    # Combine AC_info, site_info, and either AS_info or quasi_info into a single 'info' struct.
    # This is required because downstream, adjust_vcf_incompatible_types() expects all VCF INFO fields
    # to be present inside ht.info. If fields are left in separate structs (like ht.AS_info or ht.quasi_info),
    # they will not be exported to the VCF INFO column. By merging them into ht.info, we ensure all
    # relevant annotations are included in the VCF output.
    as_info_ht = build_info_ht(info_ht, 'AS_info')
    quasi_info_ht = build_info_ht(info_ht, 'quasi_info')

    logging.info(f'Writing combined AS, quasi-AS, and site-level HT to {out_ht_path}')
    info_ht.write(str(out_ht_path), overwrite=True)
    return as_info_ht, quasi_info_ht


def _filter_rows_and_add_tags(mt: hl.MatrixTable) -> hl.MatrixTable:
    # Filter to only non-reference sites.
    # An example of a variant with hl.len(mt.alleles) > 1 BUT NOT
    # hl.agg.any(mt.LGT.is_non_ref()) is a variant that spans a deletion,
    # which was however filtered out, so the LGT was set to NA, however the site
    # was preserved to account for the presence of that spanning deletion.
    # locus   alleles    LGT
    # chr1:1 ["GCT","G"] 0/1
    # chr1:3 ["T","*"]   NA
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # annotate site level DP as site_dp onto the mt rows to avoid name collision
    mt = mt.annotate_rows(site_dp=hl.agg.sum(mt.DP))

    # Add AN tag as ANS
    return mt.annotate_rows(ANS=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2)


def _create_info_ht(mt: hl.MatrixTable, n_partitions: int) -> hl.Table:
    """Create info table from vcf matrix table"""

    info_ht: hl.Table = default_compute_info(
        mt,
        as_annotations=True,
        site_annotations=True,
        n_partitions=n_partitions,
    )

    info_ht = info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))
    return info_ht
