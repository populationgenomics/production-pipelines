"""
This Stage takes a background VDS dataset (e.g. thousand genomes, HGDP, or HGDP-1KG) and densifies it to a MatrixTable
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch


def densify_background_dataset(data_path: str, qc_variants_ht: hl.Table) -> hl.MatrixTable:
    """

    Args:
        data_path ():
        qc_variants_ht ():
        tmp_path ():

    Returns:
        A densified MT
    """
    if data_path.endswith('.mt'):
        background_mt = hl.read_matrix_table(data_path)
        background_mt = hl.split_multi(background_mt, filter_changed_loci=True)
        background_mt = background_mt.semi_join_rows(qc_variants_ht)
        background_mt = background_mt.densify()
    elif data_path.endswith('.vds'):
        background_vds = hl.vds.read_vds(data_path)
        background_vds = hl.vds.split_multi(background_vds, filter_changed_loci=True)
        background_vds = hl.vds.filter_variants(background_vds, qc_variants_ht)
        background_mt = hl.vds.to_dense_mt(background_vds)
    else:
        raise ValueError('Background dataset path must be either .mt or .vds')

    return background_mt


def cli_main():
    """
    A command-line entrypoint for the densification of a background dataset
    """

    parser = ArgumentParser()
    parser.add_argument('--background-vds', help='Path to the background VDS dataset')
    parser.add_argument('--qc-variants-table', help='Path to the QC variants table')
    parser.add_argument('--dense-out', help='Path to the output dense MT')
    parser.add_argument('--tmp-path', help='Path to a temporary directory')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(
        background_vds=args.background_vds,
        qc_variants_table=args.qc_variants_table,
        dense_out=args.dense_out,
        tmp_path=args.tmp_path,
    )


def main(background_vds: str, qc_variants_table: str, dense_out: str, tmp_path: str):
    """
    Entrypoint for the densification of a background dataset
    """

    init_batch()

    background_vds = hl.vds.read_vds(background_vds)
    qc_variants_ht = hl.read_table(qc_variants_table)

    dense_background_mt = densify_background_dataset(background_vds, qc_variants_ht)

    logging.info(f'Writing dense background MT to {dense_out}')
    dense_background_mt.write(dense_out)


if __name__ == '__main__':
    cli_main()
