import logging
from typing import Optional, Union

import hail as hl

from cpg_utils import Path, to_path
from cpg_utils.config import output_path
from cpg_workflows.utils import can_reuse
from gnomad.utils import reference_genome, sparse_mt
from gnomad.utils.annotations import generate_freq_group_membership_array

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_coverage_stats(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    coverage_over_x_bins: list[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    row_key_fields: list[str] = ["locus"],
    strata_expr: Optional[list[dict[str, hl.expr.StringExpression]]] = None,
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

    :param mtds: Input sparse MT or VDS
    :param reference_ht: Input reference HT
    :param interval_ht: Optional Table containing intervals to filter to
    :param coverage_over_x_bins: List of boundaries for computing samples over X
    :param row_key_fields: List of row key fields to use for joining `mtds` with
        `reference_ht`
    :param strata_expr: Optional list of dicts containing expressions to stratify the
        coverage stats by.
    :return: Table with per-base coverage stats
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds

    if strata_expr is None:
        strata_expr = [{}]
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
    print(f'group_membership_ht: {group_membership_ht.describe()}')

    n_samples = group_membership_ht.count()
    sample_counts = group_membership_ht.index_globals().freq_meta_sample_count

    logger.info("Computing coverage stats on %d samples.", n_samples)
    # Filter datasets to interval list
    if interval_ht is not None:
        reference_ht = reference_ht.filter(
            hl.is_defined(interval_ht[reference_ht.locus]),
        )

        if is_vds:
            mtds = hl.vds.filter_intervals(
                vds=mtds,
                intervals=interval_ht,
                split_reference_blocks=True,
            )
        else:
            raise NotImplementedError(
                "Filtering to an interval list for a sparse Matrix Table is currently" " not supported.",
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

    if is_vds:
        mtds = hl.vds.VariantDataset(
            mtds.reference_data.select_entries("END", "DP").select_cols().select_rows(),
            join_with_ref(mtds.variant_data),
        )

        # Densify
        mt = hl.vds.to_dense_mt(mtds)
    else:
        mtds = join_with_ref(mtds)
        # Densify
        mt = hl.experimental.densify(mtds)

    print(f'mt: {mt.describe()}')
    # Filter rows where the reference is missing
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null
    mt = mt.unfilter_entries()

    # Annotate with group membership
    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )

    print(f'mt after annotating with group membership: {mt.describe()}')
    print(f'mt after annotating with group membership: {mt.show()}')
    print(f'group_membership_ht show(): {group_membership_ht.show()}')
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
    coverage_ht: hl.Table = compute_coverage_stats(vds, reference_ht)
    logging.info(f'coverage_ht: {coverage_ht.describe()}')

    logging.info(f'Writing coverage data to {out_path}...')
    coverage_ht.write(out_path, overwrite=True)
    logging.info('Coverage stats written to table.')
    return coverage_ht
