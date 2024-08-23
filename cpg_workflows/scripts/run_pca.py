import logging

import pandas as pd

import hail as hl

from gnomad.utils.filtering import filter_to_autosomes


def pc_project(
    mt: hl.MatrixTable,
    loadings_ht: hl.Table,
    loading_location: str = "loadings",
    af_location: str = "pca_af",
) -> hl.Table:
    """
    Project samples in `mt` on pre-computed PCs.

    :param mt: MT containing the samples to project
    :param loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param loading_location: Location of expression for loadings in `loadings_ht`
    :param af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    """
    n_variants = loadings_ht.count()

    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location],
    )

    mt = mt.filter_rows(
        hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) & (mt.pca_af > 0) & (mt.pca_af < 1),
    )

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(
        n_variants * 2 * mt.pca_af * (1 - mt.pca_af),
    )

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select("scores")


def run() -> tuple[hl.Table, hl.Table, hl.Table]:

    output_directory = 'gs://cpg-bioheart-test/pca_differences'
    hl.init(default_reference='GRCh38')
    dense_mt = hl.read_matrix_table('gs://cpg-bioheart-test/pca_differences/dense_subset_checkpoint.mt')
    relateds_to_drop_path = 'gs://cpg-bioheart-test/large_cohort/tenk10k1-0/relateds_to_drop.ht'
    related_samples_to_drop = hl.read_table(relateds_to_drop_path)
    n_pcs = 5
    autosomes_only = True

    if autosomes_only:
        logging.info('Filtering to autosomes')
        dense_mt = filter_to_autosomes(dense_mt)

    pca_mt = dense_mt

    logging.info(f'relateds_to_drop: {related_samples_to_drop.show()}')
    if related_samples_to_drop:
        pca_mt = pca_mt.filter_cols(
            hl.is_missing(related_samples_to_drop[pca_mt.col_key]),
        )
    logging.info(f'pca_mt cols after relateds_to_drop: {pca_mt.cols().show()}')
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(
        pca_mt.GT,
        k=n_pcs,
        compute_loadings=True,
    )

    pca_af_ht = pca_mt.annotate_rows(
        pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2,
    ).rows()
    pca_loadings = pca_loadings.annotate(
        pca_af=pca_af_ht[pca_loadings.key].pca_af,
    )  # TODO: Evaluate if needed to write results at this point if relateds or not

    pca_loadings = pca_loadings.persist()
    pca_scores = pca_scores.persist()
    project_pca_mt = dense_mt.filter_cols(hl.is_missing(pca_mt.cols()[dense_mt.col_key]))
    projected_scores = pc_project(project_pca_mt, pca_loadings)
    pca_scores = pca_scores.union(projected_scores)

    print("PCA evaluations, scores, and loadings computed.")

    print("Checkpointing eigenvalues...")
    eigenvalues_ht = hl.Table.from_pandas(pd.DataFrame(pca_evals, columns=['f0']))
    pca_evals_checkpoint = eigenvalues_ht.checkpoint(output_directory + '/eigenvalues2.ht', overwrite=True)

    print("Checkpointing scores...")
    pca_scores_checkpoint = pca_scores.checkpoint(output_directory + '/scores2.ht', overwrite=True)

    print("Checkpointing loadings...")
    pca_loadings_checkpoint = pca_loadings.checkpoint(output_directory + '/loadings2.ht', overwrite=True)

    return (
        pca_evals_checkpoint,
        pca_scores_checkpoint,
        pca_loadings_checkpoint,
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
