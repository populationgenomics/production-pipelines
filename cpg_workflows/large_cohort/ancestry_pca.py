import logging
import pickle
from random import sample

import pandas as pd

import hail as hl

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, reference_path
from cpg_workflows.utils import can_reuse
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds

MIN_N_PCS = 3  # for one PC1 vs PC2 plot
MIN_N_SAMPLES = 10


def add_background(
    dense_mt: hl.MatrixTable,
    sample_qc_ht: hl.Table,
    tmp_prefix: Path,
) -> tuple[hl.MatrixTable, hl.Table]:
    """
    Add background dataset samples to the dense MT and sample QC HT.
    """
    sites_table = get_config()['references']['ancestry']['sites_table']
    allow_missing_columns = get_config()['large_cohort']['pca_background'].get('allow_missing_columns', False)
    drop_columns = get_config()['large_cohort']['pca_background'].get('drop_columns')
    qc_variants_ht = hl.read_table(sites_table)
    dense_mt = dense_mt.select_cols().select_rows().select_entries('GT', 'GQ', 'DP', 'AD')
    dataset = get_config()['large_cohort']['pca_background']
    populations_to_filter = get_config()['large_cohort']['pca_background'].get('superpopulation_to_filter', False)
    for dataset in get_config()['large_cohort']['pca_background']['datasets']:
        dataset_dict = get_config()['large_cohort']['pca_background'][dataset]
        path = dataset_dict['dataset_path']
        logging.info(f'Adding background dataset {path}')
        if to_path(path).suffix == '.mt':
            background_mt = hl.read_matrix_table(str(path))
            background_mt = hl.split_multi(background_mt, filter_changed_loci=True)
            background_mt = background_mt.semi_join_rows(qc_variants_ht)
            background_mt = background_mt.densify()
        elif to_path(path).suffix == '.vds':
            background_mt_checkpoint_path = tmp_prefix / 'densified_background_mt.mt'
            if can_reuse(background_mt_checkpoint_path):
                logging.info(f'Reusing densified background mt from {background_mt_checkpoint_path}')
                background_mt = hl.read_matrix_table(str(background_mt_checkpoint_path))
            else:
                background_vds = hl.vds.read_vds(str(path))
                background_vds = hl.vds.split_multi(background_vds, filter_changed_loci=True)
                background_vds = hl.vds.filter_variants(background_vds, qc_variants_ht)
                background_mt = hl.vds.to_dense_mt(background_vds)
                logging.info(f'Checkpointing background_mt to {background_mt_checkpoint_path}')
                background_mt = background_mt.checkpoint(str(background_mt_checkpoint_path), overwrite=True)
                logging.info('Finished checkpointing densified_background_mt')
                # annotate background mt with metadata info derived from SampleQC stage
            metadata_tables = []
            for path in dataset_dict['metadata_table']:
                sample_qc_background = hl.read_table(path)
                metadata_tables.append(sample_qc_background)
            metadata_tables = hl.Table.union(*metadata_tables, unify=allow_missing_columns)
            background_mt = background_mt.annotate_cols(**metadata_tables[background_mt.col_key])
            if populations_to_filter:
                logging.info(f'Filtering background samples by {populations_to_filter}')
                background_mt = background_mt.filter_cols(
                    hl.literal(populations_to_filter).contains(background_mt.superpopulation),
                )
                logging.info(f'Finished filtering background, kept samples that are {populations_to_filter}')
        else:
            raise ValueError('Background dataset path must be either .mt or .vds')

        # save metadata info before merging dense and background datasets
        ht = background_mt.cols()
        background_mt = background_mt.select_cols().select_rows().select_entries('GT', 'GQ', 'DP', 'AD')
        background_mt = background_mt.naive_coalesce(5000)
        # combine dense dataset with background population dataset
        dense_mt = dense_mt.union_cols(background_mt)
        sample_qc_ht = sample_qc_ht.union(ht, unify=allow_missing_columns)

    if drop_columns:
        sample_qc_ht = sample_qc_ht.drop(*drop_columns)

    return dense_mt, sample_qc_ht


def run(
    dense_mt_path: Path,
    sample_qc_ht_path: Path,
    relateds_to_drop_ht_path: Path,
    tmp_prefix: Path,
    out_scores_ht_path: Path,
    out_eigenvalues_ht_path: Path,
    out_loadings_ht_path: Path,
    out_inferred_pop_ht_path: Path,
    out_sample_qc_ht_path: Path,
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
    dense_mt = hl.read_matrix_table(str(dense_mt_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))

    # If requested, subset the dense_mt and sample_qc_ht to the samples provided in the config
    sgids_remove = get_config()['large_cohort'].get('pca_samples_to_remove', [])
    if not sgids_remove:
        logging.info('No specific samples provided for removal. Continuing with the full cohort.')
    else:
        logging.info(f'Removing samples prior to PCA analysis. Removing {sgids_remove}')
        dense_mt = dense_mt.filter_cols(~hl.literal(sgids_remove).contains(dense_mt.s))
        sample_qc_ht = sample_qc_ht.filter(~hl.literal(sgids_remove).contains(sample_qc_ht.s))

    pca_background = get_config()['large_cohort'].get('pca_background', {})
    if 'datasets' in pca_background:
        dense_mt_checkpoint_path = tmp_prefix / 'modified_dense_mt.mt'
        sample_qc_ht_checkpoint_path = tmp_prefix / 'modified_sample_qc.ht'
        if can_reuse(dense_mt_checkpoint_path) and can_reuse(sample_qc_ht_checkpoint_path):
            dense_mt = hl.read_matrix_table(str(dense_mt_checkpoint_path))
            sample_qc_ht = hl.read_table(str(sample_qc_ht_checkpoint_path))
        else:
            logging.info(f'Adding background datasets using following config: {pca_background}')
            dense_mt, sample_qc_ht = add_background(dense_mt, sample_qc_ht, tmp_prefix)
            logging.info(f'Checkpointing dense_mt to {dense_mt_checkpoint_path}')
            dense_mt = dense_mt.checkpoint(str(dense_mt_checkpoint_path), overwrite=True)
            logging.info(f'Checkpointing sample_qc_ht to {sample_qc_ht_checkpoint_path}')
            sample_qc_ht = sample_qc_ht.checkpoint(str(sample_qc_ht_checkpoint_path), overwrite=True)

    logging.info(
        f'Running PCA on {dense_mt.count_cols()} samples, {dense_mt.count_rows()} sites, using {n_pcs} PCs',
    )
    scores_ht, eigenvalues_ht, loadings_ht = _run_pca_ancestry_analysis(
        mt=dense_mt,
        sample_to_drop_ht=relateds_to_drop_ht,
        n_pcs=n_pcs,
        out_scores_ht_path=out_scores_ht_path,
        out_eigenvalues_ht_path=out_eigenvalues_ht_path,
        out_loadings_ht_path=out_loadings_ht_path,
    )

    scores_ht = hl.read_table(str(out_scores_ht_path))
    training_pop_ht = sample_qc_ht.filter(hl.is_defined(sample_qc_ht['training_pop']))
    pop_ht = _infer_pop_labels(
        scores_ht=scores_ht,
        training_pop_ht=training_pop_ht,
        tmp_prefix=tmp_prefix,
        min_prob=min_pop_prob,
        n_pcs=n_pcs,
        out_ht_path=out_inferred_pop_ht_path,
    )
    sample_qc_ht.annotate(**pop_ht[sample_qc_ht.key])
    sample_qc_ht.checkpoint(str(out_sample_qc_ht_path), overwrite=True)
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
        not fp or can_reuse(fp, overwrite=True)
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

    if n_pcs < MIN_N_PCS:
        raise ValueError(f'The number of PCs must be at least {MIN_N_PCS}, got {n_pcs}')

    samples_to_use = mt.count_cols()
    logging.info(f'Total sample number in the matrix table: {samples_to_use}')
    if sample_to_drop_ht is not None:
        samples_to_drop = sample_to_drop_ht.count()
        logging.info(f'Determined {samples_to_drop} relateds to drop')
        samples_to_use -= samples_to_drop
        logging.info(
            f'Removing the {samples_to_drop} relateds from the list of samples used '
            f'for PCA, got remaining {samples_to_use} samples',
        )

    if samples_to_use < MIN_N_SAMPLES:
        raise ValueError(
            f'The number of samples after removing relateds if too low for the PCA '
            f'analysis. Got {samples_to_use}, but need at least {MIN_N_SAMPLES}',
        )

    if n_pcs > samples_to_use:
        logging.info(f'Adjusting the number of PCs not to exceed the number of samples:{n_pcs} -> {samples_to_use}')
        n_pcs = samples_to_use

    eigenvalues, scores_ht, loadings_ht = run_pca_with_relateds(mt, sample_to_drop_ht, n_pcs=n_pcs)
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
    if can_reuse(out_ht_path, overwrite=True):
        return hl.read_table(str(out_ht_path))

    if training_pop_ht.count() < 2:
        logging.warning(
            'Need at least 2 samples with known `population` label to run PCA '
            'and assign population labels to remaining samples',
        )
        pop_ht = scores_ht.annotate(
            training_pop=hl.missing(hl.tstr),
            pop='Other',
            is_training=False,
            pca_scores=scores_ht.scores,
        )
        return pop_ht.checkpoint(str(out_ht_path), overwrite=True)

    logging.info(
        'Using calculated PCA scores as well as training samples with known '
        '`population` label to assign population labels to remaining samples',
    )
    scores_ht = scores_ht.annotate(training_pop=training_pop_ht[scores_ht.key].training_pop)

    def _run_assign_population_pcs(pop_pca_scores_ht_, min_prob_):
        examples_num = pop_pca_scores_ht_.aggregate(hl.agg.count_where(hl.is_defined(pop_pca_scores_ht_.training_pop)))
        logging.info(f'Running RF using {examples_num} training examples')
        pop_ht_, pops_rf_model_ = assign_population_pcs(
            pop_pca_scores_ht_,
            pc_cols=pop_pca_scores_ht_.scores[:n_pcs],
            known_col='training_pop',
            min_prob=min_prob_,
        )
        n_mislabeled_samples_ = pop_ht_.aggregate(hl.agg.count_where(pop_ht_.training_pop != pop_ht_.pop))
        return pop_ht_, pops_rf_model_, n_mislabeled_samples_

    pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(scores_ht, min_prob)
    while n_mislabeled_samples > max_mislabeled_training_samples:
        logging.info(
            f'Found {n_mislabeled_samples} samples '
            f'labeled differently from their known pop. '
            f'Re-running without them.',
        )

        pop_ht = pop_ht[scores_ht.key]
        pop_pca_scores_ht = scores_ht.annotate(
            training_pop=hl.or_missing((pop_ht.training_pop == pop_ht.pop), scores_ht.training_pop),
        ).persist()

        pop_ht, pops_rf_model, n_mislabeled_samples = _run_assign_population_pcs(pop_pca_scores_ht, min_prob)

    # Writing a tab delimited file indicating inferred sample populations
    pop_tsv_file = tmp_prefix / 'RF_pop_assignments.txt.gz'
    if not can_reuse(pop_tsv_file, overwrite=True):
        pc_cnt = min(hl.min(10, hl.len(pop_ht.pca_scores)).collect())
        pop_ht.transmute(**{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(pc_cnt)}).export(str(pop_tsv_file))

    # Writing the RF model used for inferring sample populations
    pop_rf_file = tmp_prefix / 'pop.RF_fit.pickle'
    if not can_reuse(pop_rf_file, overwrite=True):
        with hl.hadoop_open(str(pop_rf_file), 'wb') as out:
            pickle.dump(pops_rf_model, out)

    pop_ht = pop_ht.annotate(is_training=hl.is_defined(training_pop_ht[pop_ht.key]))
    return pop_ht.checkpoint(str(out_ht_path), overwrite=True)
