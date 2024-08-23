import logging

import pandas as pd

import hail as hl

from gnomad.utils.filtering import filter_to_autosomes


def run() -> tuple[hl.Table, hl.Table, hl.Table]:

    output_directory = 'gs://cpg-bioheart-test/pca_differences'
    hl.init(default_reference='GRCh38')
    dense_mt = hl.read_matrix_table('gs://cpg-bioheart-test/pca_differences/dense_subset_checkpoint.mt')
    relateds_to_drop_path = 'gs://cpg-bioheart-test/large_cohort/tenk10k1-0/relateds_to_drop.ht'
    related_samples_to_drop = hl.read_table(relateds_to_drop_path)
    n_pcs = 5

    logging.info('Filtering to autosomes')
    dense_mt = filter_to_autosomes(dense_mt)

    logging.info(f'relateds_to_drop: {related_samples_to_drop.show()}')
    if related_samples_to_drop:
        dense_mt = dense_mt.filter_cols(
            hl.is_missing(related_samples_to_drop[dense_mt.col_key]),
        )
    logging.info(f'dense_mt cols after relateds_to_drop: {dense_mt.cols().show()}')
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(
        dense_mt.GT,
        k=n_pcs,
        compute_loadings=True,
    )

    print("PCA evaluations, scores, and loadings computed.")

    print("Checkpointing eigenvalues...")
    eigenvalues_ht = hl.Table.from_pandas(pd.DataFrame(pca_evals, columns=['f0']))
    pca_evals_checkpoint = eigenvalues_ht.checkpoint(output_directory + '/eigenvalues.ht', overwrite=True)

    print("Checkpointing scores...")
    pca_scores_checkpoint = pca_scores.checkpoint(output_directory + '/scores.ht', overwrite=True)

    print("Checkpointing loadings...")
    pca_loadings_checkpoint = pca_loadings.checkpoint(output_directory + '/loadings.ht', overwrite=True)

    return (
        pca_evals_checkpoint,
        pca_scores_checkpoint,
        pca_loadings_checkpoint,
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
