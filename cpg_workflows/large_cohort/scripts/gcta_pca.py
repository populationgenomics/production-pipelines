"""
input is dense mt
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import get_batch, init_batch


def import_plink_files(output_path: str) -> hl.MatrixTable:
    return dict(bed=f'{output_path}.bed', bim=f'{output_path}.bim', fam=f'{output_path}.fam')


def create_plink_files(dense_mt_path: str, output_path: str) -> hl.MatrixTable:
    dense_mt = hl.read_matrix_table(dense_mt_path)
    hl.export_plink(dense_mt, output_path, ind_id=dense_mt.s)
    return dict(bed=f'{output_path}.bed', bim=f'{output_path}.bim', fam=f'{output_path}.fam')


def cli_main():

    init_batch()

    parser = ArgumentParser()
    parser.add_argument('--dense-mt-path', help='Path to the dense MT')
    parser.add_argument('--output-path', help='Path to the output plink files')
    parser.add_argument('--version', help='Version of the plink files')
    parser.add_argument('--create-plink', help='Create plink files', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(
        dense_mt_path=args.dense_mt_path,
        output_path=args.output_path,
        version=args.version,
        create_plink=args.create_plink,
    )


def main(dense_mt_path: str, output_path: str, version: str, create_plink: bool | None = False):

    output_path = output_path + f'/{version}'
    if create_plink:
        plink_files = create_plink_files(dense_mt_path=dense_mt_path, output_path=output_path)
    else:
        plink_files = import_plink_files(output_path=output_path)

    b = get_batch()
    # create_plink_j = b.new_job('Create Plink files')
    # create_plink_j.image(image_path('cpg_workflows'))
    # create_plink_j.declare_resource_group(ofile={
    #     'bed': '{root}.bed',
    #     'bim': '{root}.bim',
    #     'fam': '{root}.fam',
    # })

    j = b.new_job('Create GRM')
    j.image(image_path('gcta'))
    j.command(
        f'gcta --bfile {output_path} --make-grm --out {j.out_grm}',
    )

    b.write_output(j.out_grm, f'{output_path}.grm')
    b.run()


if __name__ == '__main__':
    cli_main()
