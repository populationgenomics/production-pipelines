import logging
from re import split

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


def compute_coverage_stats(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    intervals: list[hl.Interval] = [],
    coverage_over_x_bins: list[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    row_key_fields: list[str] = ["locus"],
    strata_expr: list[dict[str, hl.expr.StringExpression]] = [],
    split_reference_blocks: bool = False,
) -> hl.Table:
    """
    Compute coverage statistics for every base of the `reference_ht` provided.

    The following coverage stats are calculated:
        - mean
        - median
        - total DP
        - fraction of samples with coverage above X, for each x in `coverage_over_x_bins`

    The `reference_ht` is a Table that contains a row for each locus coverage that should be
    computed on. It needs to be keyed by `locus`. The `reference_ht` can e.g. be
    created using `get_reference_ht`.

    :param vds: Input sparse MT or VDS
    :param reference_ht: Input reference HT
    :param intervals: Optional Table containing intervals to filter to
    :param coverage_over_x_bins: List of boundaries for computing samples over X
    :param row_key_fields: List of row key fields to use for joining `vds` with
        `reference_ht`
    :param strata_expr: Optional list of dicts containing expressions to stratify the
        coverage stats by.
    :return: Table with per-base coverage stats
    """
    mt = vds.variant_data

    if strata_expr is None:
        strata_expr = {}
        no_strata = True
    else:
        no_strata = False

    # Annotate the MT cols with each of the expressions in strata_expr and redefine
    # strata_expr based on the column HT with added annotations.
    ht = mt.annotate_cols(**{k: v for d in strata_expr for k, v in d.items()}).cols()
    strata_expr = [{k: ht[k] for k in d} for d in strata_expr]

    # Use the function for creating the frequency stratified by `freq_meta`,
    # `freq_meta_sample_count`, and `group_membership` annotations to give
    # stratification group membership info for computing coverage. By default, this
    # function returns annotations where the second element is a placeholder for the
    # "raw" frequency of all samples, where the first 2 elements are the same sample
    # set, but freq_meta startswith [{"group": "adj", "group": "raw", ...]. Use
    # `no_raw_group` to exclude the "raw" group so there is a single annotation
    # representing the full samples set. `freq_meta` is updated below to remove "group"
    # from all dicts.
    group_membership_ht = generate_freq_group_membership_array(
        ht,
        strata_expr,
        no_raw_group=True,
    )
    n_samples = group_membership_ht.count()
    sample_counts = group_membership_ht.index_globals().freq_meta_sample_count

    logger.info("Computing coverage stats on %d samples.", n_samples)
    # Filter datasets to interval list
    if intervals is not None:
        # Building reference_ht based off of interavls so don't need to filter it
        # reference_ht = reference_ht.filter(
        #     hl.is_defined(intervals[reference_ht.locus]),
        # )

        vds = hl.vds.filter_intervals(
            vds=vds,
            intervals=intervals,
            split_reference_blocks=split_reference_blocks,
        )

    # Create an outer join with the reference Table
    def join_with_ref(mt: hl.MatrixTable) -> hl.MatrixTable:
        """
        Outer join MatrixTable with reference Table.

        Add 'in_ref' annotation indicating whether a given position is found in the reference Table.

        :param mt: Input MatrixTable.
        :return: MatrixTable with 'in_ref' annotation added.
        """
        keep_entries = ["DP"]
        if "END" in mt.entry:
            keep_entries.append("END")
        if "LGT" in mt.entry:
            keep_entries.append("LGT")
        if "GT" in mt.entry:
            keep_entries.append("GT")
        mt_col_key_fields = list(mt.col_key)
        mt_row_key_fields = list(mt.row_key)
        t = mt.select_entries(*keep_entries).select_cols().select_rows()
        t = t._localize_entries("__entries", "__cols")
        t = (
            t.key_by(*row_key_fields)
            .join(
                reference_ht.key_by(*row_key_fields).select(_in_ref=True),
                how="outer",
            )
            .key_by(*mt_row_key_fields)
        )
        t = t.annotate(
            __entries=hl.or_else(
                t.__entries,
                hl.range(n_samples).map(
                    lambda x: hl.missing(t.__entries.dtype.element_type),
                ),
            ),
        )

        return t._unlocalize_entries("__entries", "__cols", mt_col_key_fields)

    keep_entries = ['END', 'DP']
    # vds objects created from Hail >= 0.2.134 have additional 'LEN' entry field
    if 'LEN' in vds.reference_data.entry:
        keep_entries.append('LEN')
    vds = hl.vds.VariantDataset(
        vds.reference_data.select_entries(*keep_entries).select_cols().select_rows(),
        join_with_ref(vds.variant_data),
    )

    # Densify
    mt = hl.vds.to_dense_mt(vds)

    # Filter rows where the reference is missing
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null
    mt = mt.unfilter_entries()

    # Annotate with group membership
    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )

    # Compute coverage stats
    coverage_over_x_bins = sorted(coverage_over_x_bins)
    max_coverage_bin = coverage_over_x_bins[-1]
    hl_coverage_over_x_bins = hl.array(coverage_over_x_bins)

    # This expression creates a counter DP -> number of samples for DP between
    # 0 and max_coverage_bin
    coverage_counter_expr = hl.agg.counter(
        hl.min(max_coverage_bin, hl.or_else(mt.DP, 0)),
    )
    mean_expr = hl.agg.mean(hl.or_else(mt.DP, 0))

    # Annotate all rows with coverage stats for each strata group.
    ht = mt.select_rows(
        coverage_stats=hl.agg.array_agg(
            lambda x: hl.agg.filter(
                x,
                hl.struct(
                    coverage_counter=coverage_counter_expr,
                    mean=hl.if_else(hl.is_nan(mean_expr), 0, mean_expr),
                    median_approx=hl.or_else(
                        hl.agg.approx_median(hl.or_else(mt.DP, 0)),
                        0,
                    ),
                    total_DP=hl.agg.sum(mt.DP),
                ),
            ),
            mt.group_membership,
        ),
    ).rows()
    ht = ht.checkpoint(hl.utils.new_temp_file("coverage_stats", "ht"))

    # This expression aggregates the DP counter in reverse order of the
    # coverage_over_x_bins and computes the cumulative sum over them.
    # It needs to be in reverse order because we want the sum over samples
    # covered by > X.
    count_array_expr = ht.coverage_stats.map(
        lambda x: hl.cumulative_sum(
            hl.array(
                # The coverage was already floored to the max_coverage_bin, so no more
                # aggregation is needed for the max bin.
                [hl.int32(x.coverage_counter.get(max_coverage_bin, 0))],
                # For each of the other bins, coverage needs to be summed between the
                # boundaries.
            ).extend(
                hl.range(hl.len(hl_coverage_over_x_bins) - 1, 0, step=-1).map(
                    lambda i: hl.sum(
                        hl.range(
                            hl_coverage_over_x_bins[i - 1],
                            hl_coverage_over_x_bins[i],
                        ).map(lambda j: hl.int32(x.coverage_counter.get(j, 0))),
                    ),
                ),
            ),
        ),
    )

    ht = ht.annotate(
        coverage_stats=hl.map(
            lambda c, g, n: c.annotate(
                **{
                    f"over_{x}": g[i] / n
                    for i, x in zip(
                        range(len(coverage_over_x_bins) - 1, -1, -1),
                        # Reverse the bin index as count_array_expr has reverse order.
                        coverage_over_x_bins,
                    )
                },
            ).drop("coverage_counter"),
            ht.coverage_stats,
            count_array_expr,
            sample_counts,
        ),
    )
    current_keys = list(ht.key)
    ht = ht.key_by(*row_key_fields).select_globals().drop(*[k for k in current_keys if k not in row_key_fields])
    if no_strata:
        # If there was no stratification, move coverage_stats annotations to the top
        # level.
        ht = ht.select(**{k: ht.coverage_stats[0][k] for k in ht.coverage_stats[0]})
    else:
        # If there was stratification, add the metadata and sample count info for the
        # stratification to the globals.
        ht = ht.annotate_globals(
            coverage_stats_meta=(
                group_membership_ht.index_globals().freq_meta.map(
                    lambda x: hl.dict(x.items().filter(lambda m: m[0] != "group")),
                )
            ),
            coverage_stats_meta_sample_count=(group_membership_ht.index_globals().freq_meta_sample_count),
        )

    return ht


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
    vds: hl.vds.VariantDataset = hl.vds.read_vds(vds_path)

    # Generate coverage table
    logging.info("Computing coverage statistics.")
    coverage_ht: hl.Table = compute_coverage_stats(
        vds,
        ref_ht_joined,
        intervals=intervals,
        split_reference_blocks=False,
    )

    return coverage_ht.checkpoint(out_path, overwrite=True)
