import logging

import click

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path, to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import Batch, get_batch, init_batch
from cpg_workflows.utils import can_reuse


def make_plink(mt: hl.MatrixTable, plink_output_path: Path):
    plink_paths = [f'{plink_output_path}.{ext}' for ext in ['bed', 'bim', 'fam']]
    if can_reuse(plink_paths):
        logging.info('Plink files already exist')
        return True
    logging.info('Exporting plink')
    return hl.export_plink(mt, str(plink_output_path), ind_id=mt.s)


def make_plink_job(dense_mt: hl.MatrixTable, plink_output_path: Path) -> hb.batch._resource.PythonResult:
    logging.info('Creating plink files')
    make_plink_job = get_batch().new_python_job('Make Plink Files', ({'tool': 'hail query'}))
    make_plink_job.image(image_path('cpg_workflows'))
    logging.info('Created plink job')
    plink_result = make_plink_job.call(make_plink, dense_mt, plink_output_path)

    return plink_result


def create_GRM(
    b: hb.Batch,
    bed_file_path: str,
    bim_file_path: str,
    fam_file_path: str,
    output_path: str,
):

    # Read in PLINK files created by Hail
    bfile = b.read_input_group(bed=bed_file_path, bim=bim_file_path, fam=fam_file_path)

    # Create GRM using GCTA
    create_GRM_j = b.new_job('Create GRM')
    create_GRM_j.image(image_path('gcta'))
    # GCTA needs chromosome names to be only numbers in .bim file
    edit_chr_cmd = (
        f"awk '{{sub(/^chr/, \"\", $1); print}}' {bfile.bim} > {bfile.bim}.tmp && mv {bfile.bim}.tmp {bfile.bim}"
    )
    create_GRM_j.command(edit_chr_cmd)
    create_GRM_j.declare_resource_group(
        ofile={
            'grm.bin': '{root}.grm.bin',
            'grm.id': '{root}.grm.id',
            'grm.N.bin': '{root}.grm.N.bin',
        },
    )
    create_GRM_j.command(
        f'gcta --bfile {bfile} --make-grm --out {create_GRM_j.ofile}',
    )

    return create_GRM_j.ofile


@click.command()
@click.option('--dense-mt-path')
@click.option('--plink-output-path')
@click.option('--grm-output-path')
@click.option('--version')
def main(dense_mt_path: Path, plink_output_path: str, grm_output_path: str, version: str):
    logging.basicConfig(level=logging.INFO)
    b = get_batch()
    init_batch()
    # Make Plink files
    plink_output_path = str(to_path(plink_output_path) / version)
    grm_output_path = str(to_path(grm_output_path) / version)
    logging.info(f'plink_output_path: {plink_output_path}')

    # logging.info('Starting batch')
    # init_batch()
    # logging.info('Batch started')
    logging.info(f'Loading dense MT from {dense_mt_path}')
    dense_mt = hl.read_matrix_table(dense_mt_path)
    logging.info('Loaded dense MT')
    plink_result = make_plink(dense_mt, plink_output_path)
    logging.info('Ran make_plink()')
    logging.info(f'plink_result: {plink_result}')

    # Create GRM
    if plink_result:
        logging.info('Creating GRM')
        create_GRM_j = create_GRM(
            b=b,
            bed_file_path=f'{plink_output_path}.bed',
            bim_file_path=f'{plink_output_path}.bim',
            fam_file_path=f'{plink_output_path}.fam',
            output_path=plink_output_path,
        )
        b.write_output(create_GRM_j.ofile, grm_output_path)
    # b = get_batch()
    logging.info('Running batch')
    b.run()


if __name__ == '__main__':
    main()
