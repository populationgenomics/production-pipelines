import logging  # noqa: I001
from datetime import datetime
from math import ceil
from typing import Optional

import hail as hl

from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_path, reference_path
from cpg_utils.hail_batch import genome_build
from cpg_workflows.utils import can_reuse, exists
from gnomad.utils.annotations import (
    generate_freq_group_membership_array,
    build_freq_stratification_list,
    qual_hist_expr,
    annotate_adj,
)

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


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


def agg_by_strata(
    mt: hl.MatrixTable,
    entry_agg_funcs,
    select_fields=None,
    group_membership_ht=None,
    entry_agg_group_membership=None,
) -> hl.Table:
    """
    Get row expression for annotations of each entry aggregation function(s) by strata.

    The entry aggregation functions are applied to the MatrixTable entries and
    aggregated. If no `group_membership_ht` (like the one returned by
    `generate_freq_group_membership_array`) is supplied, `mt` must contain a
    'group_membership' annotation that is a list of bools to aggregate the columns by.

    :param mt: Input MatrixTable.
    :param entry_agg_funcs: Dict of entry aggregation functions where the
        keys of the dict are the names of the annotations and the values are tuples
        of functions. The first function is used to transform the `mt` entries in some
        way, and the second function is used to aggregate the output from the first
        function.
    :param select_fields: Optional list of row fields from `mt` to keep on the output
        Table.
    :param group_membership_ht: Optional Table containing group membership annotations
        to stratify the aggregations by. If not provided, the 'group_membership'
        annotation is expected to be present on `mt`.
    :param entry_agg_group_membership: Optional dict indicating the subset of group
        strata in 'freq_meta' to run the entry aggregation functions on. The keys of
        the dict can be any of the keys in `entry_agg_funcs` and the values are lists
        of dicts. Each dict in the list contains the strata in 'freq_meta' to use for
        the corresponding entry aggregation function. If provided, 'freq_meta' must be
        present in `group_membership_ht` or `mt` and represent the same strata as those
        in 'group_membership'. If not provided, all entries of the 'group_membership'
        annotation will have the entry aggregation functions applied to them.
    :return: Table with annotations of stratified aggregations.
    """
    if group_membership_ht is None and 'group_membership' not in mt.col:
        raise ValueError(
            "The 'group_membership' annotation is not found in the input MatrixTable "
            "and 'group_membership_ht' is not specified.",
        )

    if select_fields is None:
        select_fields = []

    if group_membership_ht is None:
        logger.info(
            "'group_membership_ht' is not specified, using sample stratification "
            "indicated by the 'group_membership' annotation on the input MatrixTable.",
        )
        group_globals = mt.index_globals()
    else:
        logger.info(
            "'group_membership_ht' is specified, using sample stratification indicated "
            "by its 'group_membership' annotation.",
        )
        group_globals = group_membership_ht.index_globals()
        mt = mt.annotate_cols(group_membership=group_membership_ht[mt.col_key].group_membership)

    global_expr = {}
    n_groups = len(mt.group_membership.take(1)[0])
    if 'adj_groups' in group_globals:
        logger.info("Using the 'adj_groups' global annotation to determine adj filtered stratification groups.")
        global_expr['adj_groups'] = group_globals.adj_groups
    elif 'freq_meta' in group_globals:
        logger.info(
            "No 'adj_groups' global annotation found, using the 'freq_meta' global "
            'annotation to determine adj filtered stratification groups.',
        )
        global_expr['adj_groups'] = group_globals.freq_meta.map(lambda x: x.get('group', 'NA') == 'adj')
    else:
        logger.info(
            "No 'adj_groups' or 'freq_meta' global annotations found. All groups will be considered non-adj.",
        )
        global_expr['adj_groups'] = hl.range(n_groups).map(lambda x: False)

    if entry_agg_group_membership is not None and 'freq_meta' not in group_globals:
        raise ValueError(
            "The 'freq_meta' global annotation must be supplied when the 'entry_agg_group_membership' is specified.",
        )

    entry_agg_group_membership = entry_agg_group_membership or {}
    entry_agg_group_membership = {
        ann: [group_globals['freq_meta'].index(s) for s in strata] for ann, strata in entry_agg_group_membership.items()
    }

    n_adj_groups = hl.eval(hl.len(global_expr['adj_groups']))
    if n_adj_groups != n_groups:
        raise ValueError(
            f"The number of elements in the 'adj_groups' ({n_adj_groups}) global "
            'annotation does not match the number of elements in the '
            f"'group_membership' annotation ({n_groups})!",
        )

    # Keep only the entries needed for the aggregation functions.
    select_expr = {**{ann: f[0](mt) for ann, f in entry_agg_funcs.items()}}
    has_adj = hl.eval(hl.any(global_expr['adj_groups']))
    if has_adj:
        select_expr['adj'] = mt.adj

    mt = mt.select_entries(**select_expr)

    # Convert MT to HT with a row annotation that is an array of all samples entries
    # for that variant.
    ht = mt.localize_entries('entries', 'cols')

    # For each stratification group in group_membership, determine the indices of the
    # samples that belong to that group.
    global_expr['indices_by_group'] = hl.range(n_groups).map(
        lambda g_i: hl.range(mt.count_cols()).filter(lambda s_i: ht.cols[s_i].group_membership[g_i]),
    )
    ht = ht.annotate_globals(**global_expr)

    # Pull out each annotation that will be used in the array aggregation below as its
    # own ArrayExpression. This is important to prevent memory issues when performing
    # the below array aggregations.
    ht = ht.select(
        *select_fields,
        **{ann: ht.entries.map(lambda e: e[ann]) for ann in select_expr.keys()},
    )

    def _agg_by_group(
        indices_by_group_expr: hl.expr.ArrayExpression,
        adj_groups_expr: hl.expr.ArrayExpression,
        agg_func,
        ann_expr: hl.expr.ArrayExpression,
    ) -> hl.expr.ArrayExpression:
        """
        Aggregate `agg_expr` by group using the `agg_func` function.

        :param indices_by_group_expr: ArrayExpression of indices of samples in each group.
        :param adj_groups_expr: ArrayExpression indicating whether each group is adj.
        :param agg_func: Aggregation function to apply to `ann_expr`.
        :param ann_expr: Expression to aggregate by group.
        :return: Aggregated array expression.
        """
        f_no_adj = lambda i, *args: agg_func(ann_expr[i])
        if has_adj:
            f = lambda i, adj: hl.if_else(adj, hl.agg.filter(ht.adj[i], f_no_adj(i)), f_no_adj(i))
        else:
            f = f_no_adj

        return hl.map(
            lambda s_indices, adj: s_indices.aggregate(lambda i: f(i, adj)),
            indices_by_group_expr,
            adj_groups_expr,
        )

    # Add annotations for any supplied entry transform and aggregation functions.
    # Filter groups to only those in entry_agg_group_membership if specified.
    # If there are no specific entry group indices for an annotation, use ht[g]
    # to consider all groups without filtering.
    ht = ht.select(
        *select_fields,
        **{
            ann: _agg_by_group(  # type: ignore[misc]
                *[  # type: ignore[misc]
                    [ht[g][i] for i in entry_agg_group_membership.get(ann, [])] or ht[g]  # type: ignore[misc]
                    for g in ['indices_by_group', 'adj_groups']  # type: ignore[misc]
                ],  # type: ignore[misc]
                agg_func=f[1],
                ann_expr=ht[ann],
            )
            for ann, f in entry_agg_funcs.items()
        },
    )

    return ht.drop('cols')


def get_coverage_agg_func(dp_field: str = 'DP', max_cov_bin: int = 100):
    """
    Get a transformation and aggregation function for computing coverage.

    Can be used as an entry aggregation function in `compute_stats_per_ref_site`.

    :param dp_field: Depth field to use for computing coverage. Default is 'DP'.
    :param max_cov_bin: Maximum coverage bin (used when computing samples over X bin). Default is 100.
    :return: Tuple of functions to transform and aggregate coverage.
    """
    return (
        lambda t: hl.if_else(hl.is_missing(t[dp_field]) | hl.is_nan(t[dp_field]), 0, t[dp_field]),
        lambda dp: hl.struct(
            # This expression creates a counter DP -> number of samples for DP
            # between 0 and max_cov_bin.
            coverage_counter=hl.agg.counter(hl.min(max_cov_bin, dp)),
            mean=hl.bind(
                lambda mean_dp: hl.if_else(hl.is_nan(mean_dp), 0, mean_dp),
                hl.agg.mean(dp),
            ),
            median_approx=hl.or_else(hl.agg.approx_median(dp), 0),
            total_DP=hl.agg.sum(dp),
        ),
    )


def compute_coverage_stats(
    mtds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    coverage_over_x_bins=[1, 5, 10, 15, 20, 25, 30, 50, 100],
    row_key_fields=['locus'],
    strata_expr=None,
    group_membership_ht: Optional[hl.Table] = None,
    dp_field: str = 'DP',
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

    :param mtds: Input sparse MT or VDS.
    :param reference_ht: Input reference HT.
    :param interval_ht: Optional Table containing intervals to filter to.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
    :param row_key_fields: List of row key fields to use for joining `mtds` with
        `reference_ht`.
    :param strata_expr: Optional list of dicts containing expressions to stratify the
        coverage stats by. Only one of `group_membership_ht` or `strata_expr` can be
        specified.
    :param group_membership_ht: Optional Table containing group membership annotations
        to stratify the coverage stats by. Only one of `group_membership_ht` or
        `strata_expr` can be specified.
    :param dp_field: Name of sample depth field. Default is DP.
    :return: Table with per-base coverage stats.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds

    # Determine the genotype field.
    en = set(mt.entry)
    gt_field = en & {'GT'} or en & {'LGT'}
    if not gt_field:
        raise ValueError('No genotype field found in entry fields!')

    gt_field = gt_field.pop()

    # Add function to compute coverage stats.
    cov_bins = sorted(coverage_over_x_bins)
    rev_cov_bins = list(reversed(cov_bins))
    max_cov_bin = cov_bins[-1]
    cov_bins = hl.array(cov_bins)
    entry_agg_funcs = {'coverage_stats': get_coverage_agg_func(dp_field=dp_field, max_cov_bin=max_cov_bin)}

    ht = compute_stats_per_ref_site(
        mtds,
        reference_ht,
        entry_agg_funcs,
        row_key_fields=row_key_fields,
        interval_ht=interval_ht,
        entry_keep_fields=[gt_field, dp_field],
        strata_expr=strata_expr,
        group_membership_ht=group_membership_ht,
    )

    # This expression aggregates the DP counter in reverse order of the cov_bins and
    # computes the cumulative sum over them. It needs to be in reverse order because we
    # want the sum over samples covered by > X.
    def _cov_stats(cov_stat: hl.expr.StructExpression, n: hl.expr.Int32Expression) -> hl.expr.StructExpression:
        # The coverage was already floored to the max_coverage_bin, so no more
        # aggregation is needed for the max bin.
        count_expr = cov_stat.coverage_counter
        max_bin_expr = hl.int32(count_expr.get(max_cov_bin, 0))

        # For each of the other bins, coverage is summed between the boundaries.
        bin_expr = hl.range(hl.len(cov_bins) - 1, 0, step=-1)
        bin_expr = bin_expr.map(
            lambda i: hl.sum(hl.range(cov_bins[i - 1], cov_bins[i]).map(lambda j: hl.int32(count_expr.get(j, 0)))),
        )
        bin_expr = hl.cumulative_sum(hl.array([max_bin_expr]).extend(bin_expr))

        bin_expr = {f'over_{x}': bin_expr[i] / n for i, x in enumerate(rev_cov_bins)}

        return cov_stat.annotate(**bin_expr).drop('coverage_counter')

    ht_globals = ht.index_globals()
    if isinstance(ht.coverage_stats, hl.expr.ArrayExpression):
        ht = ht.select_globals(
            coverage_stats_meta=ht_globals.strata_meta.map(
                lambda x: hl.dict(x.items().filter(lambda m: m[0] != 'group')),
            ),
            coverage_stats_meta_sample_count=ht_globals.strata_sample_count,
        )
        cov_stats_expr = {
            'coverage_stats': hl.map(
                lambda c, n: _cov_stats(c, n),
                ht.coverage_stats,
                ht_globals.strata_sample_count,
            ),
        }
    else:
        cov_stats_expr = _cov_stats(ht.coverage_stats, ht_globals.sample_count)

    ht = ht.transmute(**cov_stats_expr)

    return ht


def densify_all_reference_sites(
    mtds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    row_key_fields=('locus',),
    entry_keep_fields=('GT',),
) -> hl.MatrixTable:
    """
    Densify a VariantDataset or Sparse MatrixTable at all sites in a reference Table.

    :param mtds: Input sparse MatrixTable or VariantDataset.
    :param reference_ht: Table of reference sites.
    :param interval_ht: Optional Table of intervals to filter to.
    :param row_key_fields: Fields to use as row key. Defaults to locus.
    :param entry_keep_fields: Fields to keep in entries before performing the
        densification. Defaults to GT.
    :return: Densified MatrixTable.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)

    if interval_ht is not None and not is_vds:
        raise NotImplementedError(
            'Filtering to an interval list for a sparse Matrix Table is currently not supported.',
        )

    # Filter datasets to interval list.
    if interval_ht is not None:
        reference_ht = reference_ht.filter(hl.is_defined(interval_ht[reference_ht.locus]))
        mtds = hl.vds.filter_intervals(vds=mtds, intervals=interval_ht, split_reference_blocks=False)

    entry_keep_fields = set(entry_keep_fields)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds
        entry_keep_fields.add('END')

    # Get the total number of samples.
    n_samples = mt.count_cols()
    mt_col_key_fields = list(mt.col_key)
    mt_row_key_fields = list(mt.row_key)
    ht = mt.select_entries(*entry_keep_fields).select_cols()

    # Localize entries and perform an outer join with the reference HT.
    ht = ht._localize_entries('__entries', '__cols')
    ht = ht.key_by(*row_key_fields)
    ht = ht.join(reference_ht.key_by(*row_key_fields).select(_in_ref=True), how='outer')
    ht = ht.key_by(*mt_row_key_fields)

    # Fill in missing entries with missing values for each entry field.
    ht = ht.annotate(
        __entries=hl.or_else(
            ht.__entries,
            hl.range(n_samples).map(lambda x: hl.missing(ht.__entries.dtype.element_type)),
        ),
    )

    # Unlocalize entries to turn the HT back to a MT.
    mt = ht._unlocalize_entries('__entries', '__cols', mt_col_key_fields)

    # Densify VDS/sparse MT at all sites.
    if is_vds:
        mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(mtds.reference_data.select_cols().select_rows(), mt))
    else:
        mt = hl.experimental.densify(mt)

    # Remove rows where the reference is missing.
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null.
    mt = mt.unfilter_entries()

    # Rekey by requested row key field and drop unused keys.
    mt = mt.key_rows_by(*row_key_fields)
    mt = mt.drop(*[k for k in mt_row_key_fields if k not in row_key_fields])

    return mt


def compute_stats_per_ref_site(
    mtds,
    reference_ht: hl.Table,
    entry_agg_funcs,
    row_key_fields=('locus',),
    interval_ht: Optional[hl.Table] = None,
    entry_keep_fields=None,
    row_keep_fields=None,
    entry_agg_group_membership=None,
    strata_expr=None,
    group_membership_ht=None,
    sex_karyotype_field=None,
) -> hl.Table:
    """
    Compute stats per site in a reference Table.

    :param mtds: Input sparse Matrix Table or VariantDataset.
    :param reference_ht: Table of reference sites.
    :param entry_agg_funcs: Dict of entry aggregation functions to perform on the
        VariantDataset/MatrixTable. The keys of the dict are the names of the
        annotations and the values are tuples of functions. The first function is used
        to transform the `mt` entries in some way, and the second function is used to
        aggregate the output from the first function.
    :param row_key_fields: Fields to use as row key. Defaults to locus.
    :param interval_ht: Optional table of intervals to filter to.
    :param entry_keep_fields: Fields to keep in entries before performing the
        densification in `densify_all_reference_sites`. Should include any fields
        needed for the functions in `entry_agg_funcs`. By default, only GT or LGT is
        kept.
    :param row_keep_fields: Fields to keep in rows after performing the stats
        aggregation. By default, only the row key fields are kept.
    :param entry_agg_group_membership: Optional dict indicating the subset of group
        strata in 'freq_meta' to use the entry aggregation functions on. The keys of
        the dict can be any of the keys in `entry_agg_funcs` and the values are lists
        of dicts. Each dict in the list contains the strata in 'freq_meta' to use for
        the corresponding entry aggregation function. If provided, 'freq_meta' must be
        present in `group_membership_ht` and represent the same strata as those in
        'group_membership'. If not provided, all entries of the 'group_membership'
        annotation will have the entry aggregation functions applied to them.
    :param strata_expr: Optional list of dicts of expressions to stratify by.
    :param group_membership_ht: Optional Table of group membership annotations.
    :param sex_karyotype_field: Optional field to use to adjust genotypes for sex
        karyotype before stats aggregation. If provided, the field must be present in
        the columns of `mtds` (variant_data MT if `mtds` is a VDS) and use "XX" and
        "XY" as values. If not provided, no sex karyotype adjustment is performed.
        Default is None.
    :return: Table of stats per site.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        mt = mtds

    if sex_karyotype_field is not None and sex_karyotype_field not in mt.col:
        raise ValueError(
            f"The supplied 'sex_karyotype_field', {sex_karyotype_field}, is not present in the columns of the input!",
        )

    if group_membership_ht is not None and strata_expr is not None:
        raise ValueError("Only one of 'group_membership_ht' or 'strata_expr' can be specified.")

    g = {} if group_membership_ht is None else group_membership_ht.globals
    if entry_agg_group_membership is not None and 'freq_meta' not in g:
        raise ValueError(
            "The 'freq_meta' annotation must be present in 'group_membership_ht' if "
            "'entry_agg_group_membership' is specified.",
        )

    # Determine if the adj annotation is needed. It is only needed if "adj_groups" is
    # in the globals of the group_membership_ht and any entry is True, or "freq_meta"
    # is in the globals of the group_membership_ht and any entry has "group" == "adj".
    adj = hl.eval(
        hl.any(g.get('adj_groups', hl.empty_array('bool')))
        | hl.any(g.get('freq_meta', hl.empty_array('dict<str, str>')).map(lambda x: x.get('group', 'NA') == 'adj')),
    )

    # Determine the entry fields on mt that should be densified.
    # "GT" or "LGT" is required for the genotype.
    # If the adj annotation is needed then "adj" must be present on mt, or AD/LAD, DP,
    # and GQ must be present.
    en = set(mt.entry)
    gt_field = en & {'GT'} or en & {'LGT'}
    ad_field = en & {'AD'} or en & {'LAD'}
    adj_fields = en & {'adj'} or ({'DP', 'GQ'} | ad_field) if adj else set([])

    if not gt_field:
        raise ValueError('No genotype field found in entry fields!')

    if adj and not adj_fields.issubset(en):
        raise ValueError(
            "No 'adj' found in entry fields, and one of AD/LAD, DP, and GQ is missing so adj can't be computed!",
        )

    entry_keep_fields = set(entry_keep_fields or set([])) | gt_field | adj_fields

    # Write the sex karyotype field out to a temp HT so we can annotate the field back
    # onto the MT after 'densify_all_reference_sites' removes all column annotations.
    if sex_karyotype_field is not None:
        sex_karyotype_ht = (
            mt.cols()
            .select(sex_karyotype_field)
            .checkpoint(dataset_path(suffix='compute_ref_stats/sex_karyotype.ht', category='tmp'), overwrite=True)
        )
    else:
        sex_karyotype_ht = None

    # Initialize no_strata and default strata_expr if neither group_membership_ht nor
    # strata_expr is provided.
    no_strata = group_membership_ht is None and strata_expr is None
    if no_strata:
        strata_expr = {}

    if group_membership_ht is None:
        logger.warning("'group_membership_ht' is not specified, no stats are adj filtered.")

        # Annotate the MT cols with each of the expressions in strata_expr and redefine
        # strata_expr based on the column HT with added annotations.
        ht = mt.annotate_cols(**{k: v for d in strata_expr for k, v in d.items()}).cols()
        strata_expr = [{k: ht[k] for k in d} for d in strata_expr]

        # Use 'generate_freq_group_membership_array' to create a group_membership Table
        # that gives stratification group membership info based on 'strata_expr'. The
        # returned Table has the following annotations: 'freq_meta',
        # 'freq_meta_sample_count', and 'group_membership'. By default, this
        # function returns annotations where the second element is a placeholder for the
        # "raw" frequency of all samples, where the first 2 elements are the same sample
        # set, but 'freq_meta' starts with [{"group": "adj", "group": "raw", ...]. Use
        # `no_raw_group` to exclude the "raw" group so there is a single annotation
        # representing the full samples set. Update all 'freq_meta' entries' "group"
        # to "raw" because `generate_freq_group_membership_array` will return them all
        # as "adj" since it was built for frequency computation, but for the coverage
        # computation we don't want to do any filtering.
        group_membership_ht = generate_freq_group_membership_array(ht, strata_expr, no_raw_group=True)
        group_membership_ht = group_membership_ht.annotate_globals(
            freq_meta=group_membership_ht.freq_meta.map(
                lambda x: hl.dict(x.items().map(lambda m: hl.if_else(m[0] == 'group', ('group', 'raw'), m))),
            ),
        )

    if is_vds:
        rmt = mtds.reference_data
        mtds = hl.vds.VariantDataset(
            rmt.select_entries(*((set(entry_keep_fields) & set(rmt.entry)) | {'END', 'LEN'})),
            mtds.variant_data,
        )

    mt = densify_all_reference_sites(
        mtds,
        reference_ht,
        interval_ht,
        row_key_fields,
        entry_keep_fields=entry_keep_fields,
    )

    if sex_karyotype_ht is not None:
        logger.info('Adjusting genotype ploidy based on sex karyotype.')
        gt_field = gt_field.pop()
        mt = mt.annotate_cols(sex_karyotype=sex_karyotype_ht[mt.col_key][sex_karyotype_field])
        mt = mt.annotate_entries(**{gt_field: adjusted_sex_ploidy_expr(mt.locus, mt[gt_field], mt.sex_karyotype)})

    # Annotate with adj if needed.
    if adj and 'adj' not in mt.entry:
        logger.info('Annotating the MT with adj.')
        mt = annotate_adj(mt)

    ht = agg_by_strata(
        mt,
        entry_agg_funcs,
        group_membership_ht=group_membership_ht,
        select_fields=row_keep_fields,
        entry_agg_group_membership=entry_agg_group_membership,
    )
    ht = ht.select_globals().checkpoint(
        dataset_path(suffix='compute_ref_stats/agg_stats.ht', category='tmp'),
        overwrite=True,
    )

    group_globals = group_membership_ht.index_globals()
    global_expr = {}
    if no_strata:
        # If there was no stratification, move aggregated annotations to the top
        # level.
        ht = ht.select(**{ann: ht[ann][0] for ann in entry_agg_funcs})
        global_expr['sample_count'] = group_globals.freq_meta_sample_count[0]
    else:
        # If there was stratification, add the metadata and sample count info for the
        # stratification to the globals.
        global_expr['strata_meta'] = group_globals.freq_meta
        global_expr['strata_sample_count'] = group_globals.freq_meta_sample_count

    ht = ht.annotate_globals(**global_expr)

    return ht


def get_allele_number_agg_func(gt_field: str = 'GT'):
    """
    Get a transformation and aggregation function for computing the allele number.

    Can be used as an entry aggregation function in `compute_stats_per_ref_site`.

    :param gt_field: Genotype field to use for computing the allele number.
    :return: Tuple of functions to transform and aggregate the allele number.
    """
    return lambda t: t[gt_field].ploidy, hl.agg.sum


def compute_an_and_qual_hists_per_ref_site(
    vds: hl.vds.VariantDataset,
    ref_ht: hl.Table,
    sample_qc_ht: hl.Table,
    interval_ht: Optional[hl.Table] = None,
    group_membership_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Compute allele number and quality histograms per reference site.

    :param vds: Input VDS.
    :param ref_ht: Reference HT.
    :param interval_ht: Interval HT.
    :param group_membership_ht: Group membership HT.
    :return: HT with allele number and quality histograms per reference site.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        return qual_hist_expr(
            gq_expr=qual_expr[0],
            dp_expr=qual_expr[1],
            adj_expr=qual_expr[2] == 1,
            split_adj_and_raw=True,
        )

    entry_agg_funcs = {
        'AN': get_allele_number_agg_func('LGT'),
        'qual_hists': (lambda t: [t.GQ, t.DP, t.adj], _get_hists),
    }

    logger.info('Computing allele number and histograms per reference site...')
    # Below we use just the raw group for qual hist computations because qual hists
    # has its own built-in adj filtering when adj is passed as an argument and will
    # produce both adj and raw histograms.

    vmt = vds.variant_data
    vmt = vmt.annotate_cols(sex_karyotype=sample_qc_ht[vmt.s].sex_karyotype)
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    ht = compute_stats_per_ref_site(
        vds,
        ref_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=['GQ', 'DP'],
        entry_agg_group_membership={'qual_hists': [{'group': 'raw'}]},
        sex_karyotype_field='sex_karyotype',
    )
    ht = ht.annotate(qual_hists=ht.qual_hists[0])

    return ht


def run_coverage(
    vds_path: str,
    sample_qc_ht_path: str,
    relateds_to_drop_ht_path: str,
    group_membership_out_path: str,
    coverage_out_path: str,
) -> hl.Table:
    """
    Generate coverage summary statistics for all provided sites in a VDS.
    :param vds_path: Path to the VDS.
    :param coverage_out_path: Path to save the coverage table.
    :return: Coverage Hail Table.
    """
    if can_reuse(coverage_out_path):
        logger.info(f'Reusing existing coverage table at {coverage_out_path}.')
        return None

    # Load the intervals to subset to, if provided.
    intervals_ht = None
    intervals_bed = config_retrieve(['large_cohort', 'coverage', 'intervals'], None)
    if intervals_bed is not None:
        logger.info('Preparing padded intervals.')
        intervals_ht = hl.import_bed(
            intervals_bed,
            reference_genome=genome_build(),
        )
        padding = config_retrieve(['large_cohort', 'coverage', 'interval_padding'], default=50)
        intervals_ht = adjust_interval_padding(intervals_ht, padding=padding)

    # Load the sites to compute coverage for.
    sites_ht = hl.read_table(config_retrieve(['large_cohort', 'coverage', 'sites_table']))

    # Load the VDS.
    vds: hl.vds.VariantDataset = hl.vds.read_vds(vds_path)

    # Subset to QC-pass samples only.
    sample_qc_ht = hl.read_table(sample_qc_ht_path)
    relateds_to_drop_ht = hl.read_table(relateds_to_drop_ht_path)
    samples_to_remove = (
        sample_qc_ht.filter(hl.len(sample_qc_ht.filters) > 0).select().union(relateds_to_drop_ht.select()).distinct()
    )
    vds = hl.vds.filter_samples(vds, samples_to_remove, keep=False)

    # Subset the QC table to QC-pass samples only prior to calculating group membership.
    sample_qc_ht = sample_qc_ht.filter(~hl.is_defined(samples_to_remove[sample_qc_ht.s]))

    # Prepare the group membership hail table.
    if can_reuse(group_membership_out_path):
        logger.info('Reusing group membership table')
        group_membership_ht = hl.read_table(group_membership_out_path)
    else:
        logger.info('Generating frequency group membership table.')
        group_membership_ht = generate_freq_group_membership_array(
            sample_qc_ht,
            build_freq_stratification_list(
                sex_expr=sample_qc_ht.sex_karyotype,
            ),
        )
        logger.info(f'Writing group membership hail table to {group_membership_out_path}.')
        group_membership_ht = group_membership_ht.checkpoint(group_membership_out_path, overwrite=True)

    # Compute coverage statistics.
    logger.info('Computing coverage statistics.')
    coverage_ht: hl.Table = compute_coverage_stats(
        vds,
        sites_ht,
        interval_ht=intervals_ht,
        group_membership_ht=group_membership_ht,
    )

    # Repartition.
    logger.info('Repartitioning coverage hail table.')
    coverage_ht = coverage_ht.repartition(config_retrieve(['large_cohort', 'coverage', 'n_partitions']))

    # Extract the coverage stats out to the top level.
    coverage_ht = coverage_ht.annotate(**coverage_ht.coverage_stats[0]).drop("coverage_stats")

    # Write to file.
    logger.info(f'Writing coverage hail table to {coverage_out_path}.')
    coverage_ht = coverage_ht.checkpoint(coverage_out_path, overwrite=True)

    return


def run_an_calculation(
    vds_path: str,
    sample_qc_ht_path: str,
    relateds_to_drop_ht_path: str,
    group_membership_out_path: str,
    an_out_path: str,
) -> hl.Table:
    """
    Generate allele number counts and quality histograms for all provided sites in a VDS.
    :param vds_path: Path to the VDS.
    :param sample_qc_ht_path: Path to a hail table containing population and sex_karyotype annotations.
    :param coverage_out_path: Path to save the allele numbers table.
    :return: Coverage Hail Table.
    """
    if can_reuse(an_out_path):
        logger.info(f'Reusing existing AN table at {an_out_path}.')
        return None

    # Load the intervals to subset to, if provided.
    intervals_ht = None
    intervals_bed = config_retrieve(['large_cohort', 'allele_number', 'intervals'], None)
    if intervals_bed is not None:
        logger.info('Preparing padded intervals.')
        intervals_ht = hl.import_bed(
            intervals_bed,
            reference_genome=genome_build(),
        )
        padding = config_retrieve(['large_cohort', 'allele_number', 'interval_padding'], default=50)
        intervals_ht = adjust_interval_padding(intervals_ht, padding=padding)

    # Load the sites to compute allele numbers for.
    sites_ht = hl.read_table(config_retrieve(['large_cohort', 'allele_number', 'sites_table']))

    # Load the VDS.
    vds: hl.vds.VariantDataset = hl.vds.read_vds(vds_path)

    # Subset to QC-pass samples only.
    sample_qc_ht = hl.read_table(sample_qc_ht_path)
    relateds_to_drop_ht = hl.read_table(relateds_to_drop_ht_path)
    samples_to_remove = (
        sample_qc_ht.filter(hl.len(sample_qc_ht.filters) > 0).select().union(relateds_to_drop_ht.select()).distinct()
    )
    vds = hl.vds.filter_samples(vds, samples_to_remove, keep=False)

    # Subset the QC table to QC-pass samples only prior to calculating group membership.
    sample_qc_ht = sample_qc_ht.filter(~hl.is_defined(samples_to_remove[sample_qc_ht.s]))

    # Prepare the group membership hail table.
    if can_reuse(group_membership_out_path):
        logger.info('Reusing group membership table')
        group_membership_ht = hl.read_table(group_membership_out_path)
    else:
        logger.info('Generating frequency group membership table.')
        group_membership_ht = generate_freq_group_membership_array(
            sample_qc_ht,
            build_freq_stratification_list(
                sex_expr=sample_qc_ht.sex_karyotype,
                pop_expr=sample_qc_ht.population,
            ),
        )
        logger.info(f'Writing group membership hail table to {group_membership_out_path}.')
        group_membership_ht = group_membership_ht.checkpoint(group_membership_out_path, overwrite=True)

    # Compute the AN and quality histograms per reference site in the specified intervals.
    logger.info('Generating allele number and quality histograms per provided site.')
    an_ht = compute_an_and_qual_hists_per_ref_site(
        vds,
        sites_ht,
        sample_qc_ht,
        interval_ht=intervals_ht,
        group_membership_ht=group_membership_ht,
    )

    logger.info('Repartitioning allele number hail table.')
    an_ht = an_ht.repartition(config_retrieve(['large_cohort', 'allele_number', 'n_partitions']))

    logger.info(f'Writing allele number hail table to {an_out_path}.')
    an_ht = an_ht.checkpoint(an_out_path, overwrite=True)

    return
