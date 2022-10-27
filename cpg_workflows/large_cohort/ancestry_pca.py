import logging
import pickle

import hail as hl
import pandas as pd
from cpg_utils import Path
from cpg_utils.config import get_config
from gnomad.sample_qc.ancestry import run_pca_with_relateds, assign_population_pcs

from cpg_workflows.utils import can_reuse


def run(
    dense_mt_path: Path,
    sample_qc_ht_path: Path,
    relateds_to_drop_ht_path: Path,
    tmp_prefix: Path,
    out_scores_ht_path: Path,
    out_eigenvalues_ht_path: Path,
    out_loadings_ht_path: Path,
    out_inferred_pop_ht_path: Path,
) -> tuple[hl.Table, hl.Table, hl.Table, hl.Table]:
    """
    Run PCA, and return 4 hail tables:
    * scores,
    * eigenvalues,
    * loadings,
    * original sample_ht annotated with the following:
        'training_pop': str
        'pca_scores': array<float64>
        'pop': str
        'prob_<pop>': float64 (for each population label)
    """
    min_pop_prob = get_config()['large_cohort']['min_pop_prob']
    n_pcs = get_config()['large_cohort']['n_pcs']
    mt = hl.read_matrix_table(str(dense_mt_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))

    logging.info(
        f'Running PCA on {mt.count_cols()} samples, {mt.count_rows()} sites, '
        f'using {n_pcs} PCs'
    )
    scores_ht, eigenvalues_ht, loadings_ht = _run_pca_ancestry_analysis(
        mt=mt,
        sample_to_drop_ht=relateds_to_drop_ht,
        n_pcs=n_pcs,
        out_scores_ht_path=out_scores_ht_path,
        out_eigenvalues_ht_path=out_eigenvalues_ht_path,
        out_loadings_ht_path=out_loadings_ht_path,
    )

    scores_ht = hl.read_table(str(out_scores_ht_path))
    training_pop_ht = sample_qc_ht.filter(hl.is_defined(sample_qc_ht['pop']))
    training_pop_ht = training_pop_ht.annotate(training_pop=training_pop_ht.pop)
    pop_ht = _infer_pop_labels(
        scores_ht=scores_ht,
        training_pop_ht=training_pop_ht,
        tmp_prefix=tmp_prefix,
        min_prob=min_pop_prob,
        n_pcs=n_pcs,
        out_ht_path=out_inferred_pop_ht_path,
    )
    sample_qc_ht.annotate(**pop_ht[sample_qc_ht.key])
    return scores_ht, eigenvalues_ht, loadings_ht, sample_qc_ht


def _run_pca_ancestry_analysis(
    mt: hl.MatrixTable,
    sample_to_drop_ht: hl.Table | None,
    out_scores_ht_path: Path,
    out_eigenvalues_ht_path: Path,
    out_loadings_ht_path: Path,
    n_pcs: int = 16,
) -> tuple[hl.Table, hl.Table, hl.Table]:
    """
    @param mt: variants usable for PCA analysis, combined with samples
        with known populations (HGDP, 1KG, etc)
    @param sample_to_drop_ht: table with samples to drop based on
        previous relatedness analysis. With a `rank` row field
    @param out_eigenvalues_ht_path: path to a txt file to write PCA eigenvalues
    @param out_scores_ht_path: path to write PCA scores
    @param out_loadings_ht_path: path to write PCA loadings
    @param n_pcs: maximum number of principal components
    @return:
     1. scores table with a row field:
        'scores': array<float64>
     2. eigenvalues table,
     3. loadings table.
    """
    logging.info('Running PCA ancestry analysis')
    if all(
        not fp or can_reuse(fp)
        for fp in [
            out_scores_ht_path,
            out_eigenvalues_ht_path,
            out_loadings_ht_path,
        ]
    ):
        return (
            hl.read_table(str(out_scores_ht_path)),
            hl.read_table(str(out_eigenvalues_ht_path)),
            hl.read_table(str(out_loadings_ht_path)),
        )

    # Adjusting the number of principal components not to exceed the
    # number of samples
    samples_to_drop_num = 0 if sample_to_drop_ht is None else sample_to_drop_ht.count()
    min_n_pcs = 3  # for one PC1 vs PC2 plot
    n_pcs = min(min_n_pcs, n_pcs, mt.cols().count() - samples_to_drop_num)

    logging.info(f'mt.s: {mt.s.collect()}')
    if sample_to_drop_ht:
        logging.info(f'sample_to_drop_ht.s: {sample_to_drop_ht.s.collect()}')

    eigenvalues, scores_ht, loadings_ht = run_pca_with_relateds(
        mt, sample_to_drop_ht, n_pcs=n_pcs
    )
    logging.info(f'scores_ht.s: {list(scores_ht.s.collect())}')
    logging.info(f'eigenvalues: {eigenvalues}')
    eigenvalues_ht = hl.Table.from_pandas(pd.DataFrame(eigenvalues, columns=['f0']))
    return (
        scores_ht.checkpoint(str(out_scores_ht_path), overwrite=True),
        eigenvalues_ht.checkpoint(str(out_eigenvalues_ht_path), overwrite=True),
        loadings_ht.checkpoint(str(out_loadings_ht_path), overwrite=True),
    )


def _infer_pop_labels(
    scores_ht: hl.Table,
    training_pop_ht: hl.Table,
    tmp_prefix: Path,
    out_ht_path: Path,
    min_prob: float,
    n_pcs: int,
    max_mislabeled_training_samples: int = 50,
) -> hl.Table:
    """
    Take population PCA results and training data, and run random forest
    to assign global population labels.

    @param scores_ht: output table of `_run_pca_ancestry_analysis()`
        with a row field 'scores': array<float64>
    @param training_pop_ht: table with samples with defined `training_pop`: str
    @param tmp_prefix: prefix for checkpoints and intermediate files
    @param min_prob: min probability of belonging to a given population
        for the population to be set (otherwise set to `None`)
    @param max_mislabeled_training_samples: keep rerunning until the number
        of mislabeled samples is below this number
    @param n_pcs: Number of PCs to use in the RF
    @param out_ht_path: Path to write the resulting HT table
    @return: a Table with the following row fields, including `prob_<POP>`
        probability fields for each population label:
        'training_pop': str
        'pca_scores': array<float64>
        'pop': str
        'prob_CEU': float64
        'prob_YRI': float64
        ... (prob_*: float64 for each population label)
    """
    out_ht_path = out_ht_path or tmp_prefix / 'inferred_pop.ht'
    if can_reuse(out_ht_path):
        return hl.read_table(str(out_ht_path))

    if training_pop_ht.count() < 2:
        logging.warning(
            'Need at least 2 samples with known `population` label to run PCA'
            'and assign population labels to remaining samples'
        )
        pop_ht = scores_ht.annotate(
            pop='Other',
            is_training=False,
            pca_scores=scores_ht.scores,
        )
        return pop_ht.checkpoint(str(out_ht_path), overwrite=True)

    logging.info(
        'Using calculated PCA scores as well as training samples with known '
        '`population` label to assign population labels to remaining samples'
    )
    scores_ht = scores_ht.annotate(
        training_pop=training_pop_ht[scores_ht.key].training_pop
    )

    def _run_assign_population_pcs(pop_pca_scores_ht_, min_prob_):
        examples_num = pop_pca_scores_ht_.aggregate(
            hl.agg.count_where(hl.is_defined(pop_pca_scores_ht_.training_pop))
        )
        logging.info(f'Running RF using {examples_num} training examples')
        pop_ht_, pops_rf_model_ = assign_population_pcs(
            pop_pca_scores_ht_,
            pc_cols=pop_pca_scores_ht_.scores[:n_pcs],
            known_col='training_pop',
            min_prob=min_prob_,
        )
        n_mislabeled_samples_ = pop_ht.aggregate(
            hl.agg.count_where(pop_ht.training_pop != pop_ht.pop)
        )
        return pop_ht_, pops_rf_model_, n_mislabeled_samples_

    pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(
        scores_ht, min_prob
    )
    while n_mislabeled_samples > max_mislabeled_training_samples:
        logging.info(
            f'Found {n_mislabeled_samples} samples '
            f'labeled differently from their known pop. '
            f'Re-running without them.'
        )

        pop_ht = pop_ht[scores_ht.key]
        pop_pca_scores_ht = scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht.pop), scores_ht.training_pop
            )
        ).persist()

        pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(
            pop_pca_scores_ht, min_prob
        )

    # Writing a tab delimited file indicating inferred sample populations
    pop_tsv_file = tmp_prefix / 'RF_pop_assignments.txt.gz'
    if not can_reuse(pop_tsv_file):
        pc_cnt = min(hl.min(10, hl.len(pop_ht.pca_scores)).collect())
        pop_ht.transmute(
            **{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(pc_cnt)}
        ).export(pop_tsv_file)

    # Writing the RF model used for inferring sample populations
    pop_rf_file = tmp_prefix / 'pop.RF_fit.pickle'
    if not can_reuse(pop_rf_file):
        with hl.hadoop_open(pop_rf_file, 'wb') as out:
            pickle.dump(pops_rf_model, out)

    pop_ht = pop_ht.annotate(is_training=hl.is_defined(training_pop_ht[pop_ht.key]))
    return pop_ht.checkpoint(str(out_ht_path), overwrite=True)
