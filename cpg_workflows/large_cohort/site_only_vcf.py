#!/usr/bin/env python

"""
Convert to site-only table and annotate with AS fields.
"""
import logging

import hail as hl
from cpg_utils import Path
from cpg_workflows.utils import can_reuse
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.utils.sparse_mt import default_compute_info


def run(
    vds_path: Path,
    sample_qc_ht_path: Path,
    relateds_to_drop_ht_path: Path,
    out_vcf_path: Path,
    tmp_prefix: Path,
):
    vds = hl.vds.read_vds(str(vds_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))

    site_only_ht_path = tmp_prefix / 'site_only.ht'
    site_only_ht = vds_to_site_only_ht(
        vds=vds,
        sample_qc_ht=sample_qc_ht,
        relateds_to_drop_ht=relateds_to_drop_ht,
        out_ht_path=site_only_ht_path,
    )
    logging.info(f'Writing site-only VCF to {out_vcf_path}')
    assert out_vcf_path.suffix == '.bgz'
    hl.export_vcf(site_only_ht, str(out_vcf_path), tabix=True)


def vds_to_site_only_ht(
    vds: hl.vds.VariantDataset,
    sample_qc_ht: hl.Table,
    relateds_to_drop_ht: hl.Table,
    out_ht_path: Path,
) -> hl.Table:
    """
    Convert VDS into sites-only VCF-ready table.
    """
    if can_reuse(out_ht_path):
        return hl.read_table(str(out_ht_path))

    mt = hl.vds.to_dense_mt(vds)
    mt = mt.filter_cols(hl.len(sample_qc_ht[mt.col_key].filters) > 0, keep=False)
    mt = mt.filter_cols(hl.is_defined(relateds_to_drop_ht[mt.col_key]), keep=False)
    mt = _filter_rows_and_add_tags(mt)
    var_ht = _create_info_ht(mt, n_partitions=mt.n_partitions())
    var_ht = adjust_vcf_incompatible_types(
        var_ht,
        # with default INFO_VCF_AS_PIPE_DELIMITED_FIELDS, AS_VarDP will be converted
        # into a pipe-delimited value e.g.: VarDP=|132.1|140.2
        # which breaks VQSR parser (it doesn't recognise the delimiter and treats
        # it as an array with a single string value "|132.1|140.2", leading to
        # an IndexOutOfBound exception when trying to access value for second allele)
        pipe_delimited_annotations=[],
    )
    var_ht.write(str(out_ht_path), overwrite=True)
    return var_ht


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
    info_ht = default_compute_info(mt, site_annotations=True, n_partitions=n_partitions)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp)
    )
    return info_ht
