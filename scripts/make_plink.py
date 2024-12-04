import logging

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path, to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import Batch, get_batch


def export_plink(dense_mt: hl.MatrixTable, output_path: str):
    return hl.export_plink(dense_mt, output_path, ind_id=dense_mt.s)


def make_plink_job(dense_mt: hl.MatrixTable, plink_output_path: Path) -> hb.batch._resource.PythonResult:
    logging.info('Creating plink files')
    make_plink_job = get_batch().new_python_job('Make Plink Files', ({'tool': 'hail query'}))
    make_plink_job.image(image_path('cpg_workflows'))
    logging.info('Created plink job')
    plink_result = make_plink_job.call(export_plink, dense_mt, plink_output_path)

    return plink_result
