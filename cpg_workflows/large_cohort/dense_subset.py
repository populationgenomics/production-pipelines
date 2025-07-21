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

    logging.info('Filtering variants to predetermined QC variants...')
    sites_table = get_config()['references']['ancestry']['sites_table']
    qc_variants_ht = hl.read_table(sites_table)
    vds = hl.vds.filter_variants(vds, qc_variants_ht)

    if 'GT' not in vds.variant_data.entry:
        logging.info('Converting LGT to GT annotations...')
        vds.variant_data = vds.variant_data.annotate_entries(
            GT=hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA),
        )

    if 'AD' not in vds.variant_data.entry:
        logging.info('Converting LAD to AD annotations...')
        vds.variant_data = vds.variant_data.annotate_entries(
            AD=hl.vds.local_to_global(
                vds.variant_data.LAD,
                vds.variant_data.LA,
                n_alleles=hl.len(vds.variant_data.alleles),
                fill_value=0,
                number='R',
            ),
        )

    logging.info('Densifying data...')
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.select_entries('GT', 'GQ', 'DP', 'AD')

    n_partitions = get_config()['large_cohort']['dense_subset_partitions']
    mt = mt.repartition(n_partitions)
    mt = mt.checkpoint(out_dense_mt_path, overwrite=True)
    logging.info(f'Number of predetermined QC variants found in the VDS: {mt.count_rows()}')
    return mt
