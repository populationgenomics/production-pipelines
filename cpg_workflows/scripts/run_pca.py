import hail as hl


def run(
    dense_mt: hl.MatrixTable,
    output_directory: str,
    n_pcs: int = 10,
) -> tuple[hl.Table, hl.Table, hl.Table]:

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
