import logging
from typing import Optional, Union

import hail as hl

from cpg_utils import Path, to_path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.utils import can_reuse
from gnomad.utils import reference_genome, sparse_mt
from gnomad.utils.annotations import generate_freq_group_membership_array

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def merge_coverage_tables(
    coverage_table_paths: list[str],
    out_path: str,
) -> hl.Table:
    """
    Merge coverage tables.

    :param coverage_tables: List of coverage tables.
    :return: Merged coverage table.
    """
    coverage_tables = [hl.read_table(coverage_table_path) for coverage_table_path in coverage_table_paths]
    merged_coverage_table = hl.Table.union(*coverage_tables)
    return merged_coverage_table.checkpoint(out_path, overwrite=True)


def shard_vds(vds_path: str, output_dict: dict[str, str]) -> hl.vds.VariantDataset:
    """
    Shard a VDS file.
    """
    vds = hl.vds.read_vds(vds_path)
    for contig, out_path in output_dict.items():
        sharded_vds = hl.vds.filter_intervals(vds, [hl.parse_locus_interval(contig, reference_genome='GRCh38')])
        sharded_vds.write(out_path, overwrite=True)
    return sharded_vds


def run(vds_path: str, reference_ht_path: str, out_path: str) -> hl.Table:
    vds: hl.vds.VariantDataset = hl.vds.read_vds(vds_path)
    reference_ht: hl.Table = hl.read_table(reference_ht_path)

    coverage_ht: hl.Table = sparse_mt.compute_coverage_stats(vds, reference_ht)

    return coverage_ht.checkpoint(out_path, overwrite=True)


def generate_reference_coverage_ht(ref: str, contig: str, out_path: str) -> hl.ReferenceGenome:
    from cpg_utils.hail_batch import init_batch

    init_batch()
    rg = hl.get_reference(ref)
    ht: hl.Table = reference_genome.get_reference_ht(rg, [contig])
    return ht.checkpoint(out_path, overwrite=True)
