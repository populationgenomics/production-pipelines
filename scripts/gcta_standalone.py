import click

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_plink(mt: hl.MatrixTable, plink_output_path: Path):
    return hl.export_plink(mt, plink_output_path, ind_id=mt.s)


def make_plink_job(dense_mt: hl.MatrixTable, plink_output_path: Path) -> hb.batch._resource.PythonResult:
    make_plink_job = get_batch().new_python_job('Make Plink Files', ({'tool': 'hail query'}))
    make_plink_job.image(image_path('cpg_workflows'))
    plink_result = make_plink_job.call(make_plink, dense_mt, plink_output_path)

    return plink_result


@click.command()
@click.option('--dense-mt-path')
@click.option('--plink-output-path')
def main(dense_mt_path: Path, plink_output_path: Path):

    init_batch()
    dense_mt = hl.read_matrix_table(dense_mt_path)
    plink_result = make_plink_job(dense_mt, plink_output_path)

    b = get_batch()
    b.write_output(plink_result, plink_output_path)


if __name__ == '__main__':
    main()
