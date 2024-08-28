"""
This Stage is a thin wrapper around the gnomad methods PCA call
"""

import logging
from argparse import ArgumentParser

import pandas as pd

import hail as hl

from cpg_utils.config import config_retrieve
from gnomad.sample_qc.ancestry import run_pca_with_relateds

MIN_N_PCS = 3  # for one PC1 vs PC2 plot
MIN_N_SAMPLES = 10


def cli_main():
    """
    A command-line entrypoint for the ancestry add-background process
    """

    parser = ArgumentParser()
    parser.add_argument('--dense_mt', help='The Dense MT to read in')
    parser.add_argument('--related', help='The HT from Relatedness ')
    parser.add_argument('--scores_out', help='scores table output path')
    parser.add_argument('--eigen_out', help='eigenvalues table output path')
    parser.add_argument('--loadings_out', help='loadings table output path')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(
        dense_mt=args.dense_mt,
        related=args.related,
        scores_out=args.scores_out,
        eigen_out=args.eigen_out,
        loadings_out=args.loadings_out,
    )


def main(dense_mt: str, related: str, scores_out: str, eigen_out: str, loadings_out: str):
    """

    Args:
        dense_mt ():
        related ():
        scores_out (str): local path to write the scores table to
        eigen_out (): local path to write the eigenvalues table to
        loadings_out (): local path to write the loadings table to
    """

    ncpu = config_retrieve(['Ancestry', 'cores'], 8)
    hl.context.init_spark(master=f'local[{ncpu}]')
    hl.default_reference('GRCh38')
    logging.info(f'Hail version: {hl.version()}')

    mt = hl.read_matrix_table(dense_mt)
    sample_to_drop_ht = hl.read_table(related)

    n_pcs = config_retrieve(['large_cohort', 'n_pcs'], 16)

    if n_pcs < MIN_N_PCS:
        raise ValueError(f'The number of PCs must be at least {MIN_N_PCS}, got {n_pcs}')

    samples_to_use = mt.count_cols()
    logging.info(f'Total sample number in the matrix table: {samples_to_use}')

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
        n_pcs = 8
        logging.info(f'Using {n_pcs} PCs for PCA')

    eigenvalues, scores_ht, loadings_ht = run_pca_with_relateds(mt, sample_to_drop_ht, n_pcs=n_pcs)
    logging.info(f'scores_ht.s: {list(scores_ht.s.collect())}')
    logging.info(f'eigenvalues: {eigenvalues}')
    eigenvalues_ht = hl.Table.from_pandas(pd.DataFrame(eigenvalues, columns=['f0']))
    scores_ht.write(scores_out)
    eigenvalues_ht.write(eigen_out)
    loadings_ht.write(loadings_out)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()
