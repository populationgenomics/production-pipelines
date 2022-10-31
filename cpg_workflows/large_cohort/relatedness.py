import logging
import hail as hl
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop


def run(
    dense_mt_path: Path,
    sample_qc_ht_path: Path,
    out_relatedness_ht_path: Path,
    out_relateds_to_drop_ht_path: Path,
    tmp_prefix: Path,
):
    dense_mt = hl.read_matrix_table(str(dense_mt_path))
    relatedness_ht = pcrelate(
        dense_mt=dense_mt,
        out_relatedness_ht_path=out_relatedness_ht_path,
        tmp_prefix=tmp_prefix,
    )
    flag_related(
        relatedness_ht=relatedness_ht,
        sample_qc_ht=hl.read_table(str(sample_qc_ht_path)),
        out_relateds_to_drop_ht_path=out_relateds_to_drop_ht_path,
        tmp_prefix=tmp_prefix,
    )


def pcrelate(
    dense_mt: hl.MatrixTable,
    out_relatedness_ht_path: Path,
    tmp_prefix: Path,
) -> hl.Table:
    """
    Writes table with the following structure:
    Row fields:
        'i': str
        'j': str
        'kin': float64
        'ibd0': float64
        'ibd1': float64
        'ibd2': float64
    Key: ['i', 'j']
    """
    mt = dense_mt

    if can_reuse(out_relatedness_ht_path):
        return hl.read_table(str(out_relatedness_ht_path))

    mt = mt.select_entries('GT')

    logging.info('Running relatedness check')
    scores_ht_path = tmp_prefix / 'pcrelate' / 'relatedness_pca_scores.ht'
    if can_reuse(scores_ht_path):
        scores_ht = hl.read_table(str(scores_ht_path))
    else:
        sample_num = mt.cols().count()
        _, scores_ht, _ = hl.hwe_normalized_pca(
            mt.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False
        )
        scores_ht.checkpoint(str(scores_ht_path), overwrite=True)

    relatedness_ht = hl.pc_relate(
        mt.GT,
        min_individual_maf=0.01,
        scores_expr=scores_ht[mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
    )

    # Converting keys for type struct{str} to str to align
    # with the rank_ht `s` key:
    relatedness_ht = relatedness_ht.key_by(i=relatedness_ht.i.s, j=relatedness_ht.j.s)
    return relatedness_ht.checkpoint(str(out_relatedness_ht_path), overwrite=True)


def flag_related(
    relatedness_ht: hl.Table,
    sample_qc_ht: hl.Table,
    out_relateds_to_drop_ht_path: Path,
    tmp_prefix: Path,
) -> hl.Table:
    """
    Rank samples and flag samples to drop so there is only one sample per family
    left, with the highest rank in the family.

    `sample_qc_ht` has to have a `filters` and `autosomal_mean_dp` columns.
    """
    logging.info(f'Flagging related samples to drop')
    if can_reuse(out_relateds_to_drop_ht_path):
        return hl.read_table(str(out_relateds_to_drop_ht_path))

    rankings_ht_path = tmp_prefix / 'relatedness' / f'samples_rankings.ht'
    if can_reuse(rankings_ht_path):
        rank_ht = hl.read_table(str(rankings_ht_path))
    else:
        rank_ht = _compute_sample_rankings(
            ht=sample_qc_ht,
        ).checkpoint(str(rankings_ht_path), overwrite=True)

    try:
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect(rank_ht.s))
            )
        )
    except hl.ExpressionException:
        # Hail doesn't handle it with `aggregate` when none of
        # the samples is 'filtered'
        filtered_samples = hl.empty_array(hl.tstr)

    to_drop_ht = compute_related_samples_to_drop(
        relatedness_ht,
        rank_ht,
        kin_threshold=get_config()['large_cohort']['max_kin'],
        filtered_samples=filtered_samples,
    )
    to_drop_ht = to_drop_ht.checkpoint(
        str(out_relateds_to_drop_ht_path), overwrite=True
    )
    return to_drop_ht


def _compute_sample_rankings(ht: hl.Table) -> hl.Table:
    """
    Orders samples by hard filters and coverage and adds rank, which is the lower,
    the better.

    @param ht: table with a `autosomal_mean_dp` and `filters` fields.
    @return: table ordered by rank, with the following row fields:
        `rank`, `filtered`
    """
    ht = ht.drop(*list(ht.globals.dtype.keys()))
    ht = ht.select(
        'autosomal_mean_dp',
        filtered=hl.len(ht.filters) > 0,
    )
    ht = ht.order_by(ht.filtered, hl.desc(ht.autosomal_mean_dp))
    ht = ht.add_index(name='rank')
    return ht.key_by('s').select('filtered', 'rank')
