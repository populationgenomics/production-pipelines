"""
input is dense mt
"""

import logging
import struct
from argparse import ArgumentParser

import numpy as np

import hail as hl
import hailtop.batch as hb

from cpg_utils.config import image_path
from cpg_utils.hail_batch import command, get_batch, init_batch
from cpg_workflows.utils import can_reuse


def cli_main():

    init_batch()

    parser = ArgumentParser()
    parser.add_argument('--dense-mt-path', help='Path to the dense MT')
    parser.add_argument('--output-path', help='Path to the output plink files')
    parser.add_argument('--version', help='Version of the plink files')
    parser.add_argument('--n-pcs', help='Number of PCs to compute', default=10)
    parser.add_argument('--relateds-to-drop', help='HT of relateds to drop', default=None)
    parser.add_argument('--create-plink', help='Create plink files', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    create_GRM(
        dense_mt_path=args.dense_mt_path,
        output_path=args.output_path,
        version=args.version,
        n_pcs=args.n_pcs,
        relateds_to_drop=args.relateds_to_drop,
        create_plink=args.create_plink,
    )


def create_GRM(
    b: hb.Batch,
    output_path: str,
):

    # Read in PLINK files created by Hail
    logging.info(f'output_path: {output_path}')
    bfile = b.read_input_group(bed=f'{output_path}.bed', bim=f'{output_path}.bim', fam=f'{output_path}.fam')

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

    b.write_output(create_GRM_j.ofile, f'{output_path}')

    return create_GRM_j


if __name__ == '__main__':
    cli_main()
