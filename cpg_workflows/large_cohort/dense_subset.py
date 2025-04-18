import logging
from pathlib import Path

import hail as hl

from cpg_utils.config import get_config, reference_path
from cpg_workflows.utils import can_reuse


def run(vds_path: str, out_dense_mt_path: str) -> hl.MatrixTable:
    """
    Filter a sparse VariantDataset to a set of predetermined QC sites
    and return a dense MatrixTable with split multiallelics.
    @return: filtered and densified MatrixTable.
    """
    if can_reuse(out_dense_mt_path, overwrite=True):
        return hl.read_matrix_table(out_dense_mt_path)

    vds = hl.vds.read_vds(vds_path)

    logging.info('Splitting multiallelics')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logging.info('Filtering variants to predetermined QC variants...')
    sites_table = get_config()['references']['ancestry']['sites_table']
    qc_variants_ht = hl.read_table(sites_table)
    vds = hl.vds.filter_variants(vds, qc_variants_ht)
    logging.info('Densifying data...')
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.select_entries('GT', 'GQ', 'DP', 'AD')
    mt = mt.naive_coalesce(5000)
    mt = mt.checkpoint(out_dense_mt_path, overwrite=True)
    logging.info(f'Number of predetermined QC variants found in the VDS: {mt.count_rows()}')
    return mt
