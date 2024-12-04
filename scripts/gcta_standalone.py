import logging

import click

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_plink(mt: hl.MatrixTable, plink_output_path: Path):
    logging.info('Exporting plink')
    return hl.export_plink(mt, plink_output_path, ind_id=mt.s)


def make_plink_job(dense_mt: hl.MatrixTable, plink_output_path: Path) -> hb.batch._resource.PythonResult:
    logging.info('Creating plink files')
    make_plink_job = get_batch().new_python_job('Make Plink Files', ({'tool': 'hail query'}))
    make_plink_job.image(image_path('cpg_workflows'))
    logging.info('Created plink job')
    plink_result = make_plink_job.call(make_plink, dense_mt, plink_output_path)

    return plink_result


@click.command()
@click.option('--dense-mt-path')
@click.option('--plink-output-path')
def main(dense_mt_path: Path, plink_output_path: Path):
    logging.basicConfig(level=logging.INFO)

    logging.info('Starting batch')
    init_batch()
    logging.info('Batch started')
    logging.info(f'Loading dense MT from {dense_mt_path}')
    dense_mt = hl.read_matrix_table(dense_mt_path)
    logging.info('Loaded dense MT')
    plink_result = make_plink(dense_mt, plink_output_path)
    logging.info('Ran make_plink()')
    logging.info(f'plink_result: {plink_result}')

    b = get_batch()
    logging.info('Writing output')
    b.write_output(plink_result, plink_output_path)


if __name__ == '__main__':
    main()
