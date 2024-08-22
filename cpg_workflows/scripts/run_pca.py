import hail as hl


def run() -> tuple[hl.Table, hl.Table, hl.Table]:

    output_directory = 'gs://cpg-bioheart-test/pca_differences'
    hl.init(default_reference='GRCh38')
    dense_mt = hl.read_matrix_table('gs://cpg-bioheart-test/large_cohort/tenk10k1-0/dense_subset.mt')
    dense_mt_subset = dense_mt.head(1000)
    n_pcs = dense_mt_subset.count_cols()

    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(
        dense_mt.GT,
        k=n_pcs,
        compute_loadings=True,
    )
    return (
        pca_evals.checkpoint(output_directory + '/eigenvalues.ht', overwrite=True),
        pca_scores.checkpoint(output_directory + '/scores.ht', overwrite=True),
        pca_loadings.checkpoint(output_directory + '/loadings.ht', overwrite=True),
    )


if __name__ == '__main__':
    run()
