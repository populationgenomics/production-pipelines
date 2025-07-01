import logging

import hail as hl
import hailtop.batch as hb

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import genome_build
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


def adjust_interval_padding(ht: hl.Table, padding: int) -> hl.Table:
    """
    Function copied from `gnomad_qc` v4
    https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v4/annotations/compute_coverage.py#L153

    Adjust interval padding in HT.

    .. warning::

        This function can lead to overlapping intervals, so it is not recommended for
        most applications. For example, it can be used to filter a variant list to all
        variants within the returned interval list, but would not work for getting an
        aggregate statistic for each interval if the desired output is independent
        statistics.

    :param ht: HT to adjust.
    :param padding: Padding to use.
    :return: HT with adjusted interval padding.
    """
    return ht.key_by(
        interval=hl.locus_interval(
            ht.interval.start.contig,
            ht.interval.start.position - padding,
            ht.interval.end.position + padding,
            # Include the end of the intervals to capture all variants.
            reference_genome=ht.interval.start.dtype.reference_genome,
            includes_end=True,
        ),
    )


def run(
    vds_path: str,
    interval_list: list[hb.ResourceFile],
    out_path: str,
) -> hl.Table:
    """
    Generate a coverage table for a given VDS and interval.
    :param vds_path: Path to the VDS.
    :interval_list: List of .interval_list files in the form of Hail Batch ResourceFiles objects.
    :param out_path: Path to save the coverage table.
    :return: Coverage Hail Table.
    """
    from cpg_utils.hail_batch import init_batch
    from gnomad.utils.reference_genome import add_reference_sequence

    init_batch()
    rg: hl.ReferenceGenome = hl.get_reference(genome_build())
    add_reference_sequence(rg)

    # Generate reference coverage table
    intervals_ht = hl.import_locus_intervals(interval_list, reference_genome=rg)

    if config_retrieve(['workflow', 'sequencing_type']) == 'exome':
        logging.info('Adjusting interval padding for exome sequencing.')
        padding = config_retrieve(['workflow', 'exome_interval_padding'], default=50)
        intervals_ht = adjust_interval_padding(intervals_ht, padding=padding)

    # .interval_list files are imported as a Table with two columns: 'interval' and 'target'
    intervals = intervals_ht.interval.collect()
    ref_tables = []

    for interval in intervals:
        logging.info(f"Generating reference coverage table for interval: {interval}")
        ref_ht = hl.utils.range_table(
            (interval.end.position - interval.start.position),
        )
        locus_expr = hl.locus(
            contig=interval.start.contig,
            pos=ref_ht.idx + interval.start.position,
            reference_genome=rg,
        )
        ref_allele_expr = locus_expr.sequence_context().lower()
        alleles_expr = [ref_allele_expr]
        ref_ht = ref_ht.select(locus=locus_expr, alleles=alleles_expr).key_by("locus", "alleles").drop("idx")
        ref_ht = ref_ht.filter(ref_ht.alleles[0] == "n", keep=False)
        ref_tables.append(ref_ht)

    ref_ht_joined = hl.Table.union(*ref_tables)

    logging.info(f"Reading VDS from {vds_path} with intervals: {intervals}")
    vds: hl.vds.VariantDataset = hl.vds.read_vds(
        vds_path,
        intervals=intervals,
    )

    # Generate coverage table
    logging.info("Computing coverage statistics.")
    coverage_ht: hl.Table = sparse_mt.compute_coverage_stats(vds, ref_ht_joined)

    return coverage_ht.checkpoint(out_path, overwrite=True)
