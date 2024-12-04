import logging

import click
from make_plink import export_plink

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path, to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import Batch, get_batch, init_batch


def make_plink(mt: hl.MatrixTable, plink_output_path: Path):
    logging.info('Exporting plink')
    return hl.export_plink(mt, plink_output_path, ind_id=mt.s)


@click.command()
@click.option('--dense-mt-path')
@click.option('--plink-output-path')
@click.option('--version')
def main(dense_mt_path: Path, plink_output_path: str, version: str):
    logging.basicConfig(level=logging.INFO)
    b = get_batch()
    init_batch()
    # Make Plink files
    plink_output_path = to_path(plink_output_path) / version
    logging.info(f'plink_output_path: {plink_output_path}')

    # logging.info('Starting batch')
    # init_batch()
    # logging.info('Batch started')
    logging.info(f'Loading dense MT from {dense_mt_path}')
    dense_mt = hl.read_matrix_table(dense_mt_path)
    logging.info('Loaded dense MT')
    plink_result = export_plink(dense_mt, plink_output_path)
    logging.info('Ran make_plink()')
    logging.info(f'plink_result: {plink_result}')

    # Create GRM

    # b = get_batch()
    logging.info('Writing output')
    b.write_output(plink_result, str(plink_output_path))
    b.run()


if __name__ == '__main__':
    main()
