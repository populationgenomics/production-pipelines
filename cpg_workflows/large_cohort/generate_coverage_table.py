import logging

import hail as hl

from cpg_utils import Path, to_path
from cpg_utils.config import output_path
from cpg_workflows.utils import can_reuse
from gnomad.utils import reference_genome, sparse_mt


def get_reference_genome(ref_genome: str) -> hl.ReferenceGenome:
    rg = hl.get_reference(ref_genome)
    reference_ht = reference_genome.get_reference_ht(rg)
    return reference_ht


def calculate_coverage_ht(vds_path: str, out_path: str, tmp_prefix: str) -> hl.Table:
    """
    Calculate coverage for each sample.
    """
    # The `reference_ht` is a Table that contains a row for each locus coverage that should be
    # computed on. It needs to be keyed by `locus`.
    vds = hl.vds.read_vds(vds_path)

    logging.info('Calculating coverage stats...')
    reference_ht: hl.Table = get_reference_genome('GRCh38')
    # if can_reuse(str(to_path(tmp_prefix) / 'reference.ht')):
    #     logging.info(f'Reading reference_ht from {str(to_path(tmp_prefix) / "reference.ht")}...')
    #     reference_ht = hl.read_table(output_path('reference.ht', 'tmp'))
    # else:
    #     logging.info(f'Checkpointing reference_ht to {str(to_path(tmp_prefix) / "reference.ht")}...')
    #     reference_ht = reference_ht.checkpoint(str(to_path(tmp_prefix) / 'reference.ht'), overwrite=True)
    logging.info(f'reference_ht: {reference_ht.describe()}')
    coverage_ht: hl.Table = sparse_mt.compute_coverage_stats(vds, reference_ht)
    logging.info(f'coverage_ht: {coverage_ht.describe()}')

    logging.info(f'Writing coverage data to {out_path}...')
    coverage_ht.write(out_path, overwrite=True)
    logging.info('Coverage stats written to table.')
    return coverage_ht
