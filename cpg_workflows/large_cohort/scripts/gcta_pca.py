"""
input is dense mt
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import get_batch, init_batch


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
    main(
        dense_mt_path=args.dense_mt_path,
        output_path=args.output_path,
        version=args.version,
        n_pcs=args.n_pcs,
        relateds_to_drop=args.relateds_to_drop,
        create_plink=args.create_plink,
    )


def main(
    dense_mt_path: str,
    output_path: str,
    version: str,
    n_pcs: int,
    relateds_to_drop: str,
    create_plink: bool | None = False,
):

    output_path = output_path + f'/{version}'
    dense_mt = hl.read_matrix_table(dense_mt_path)
    if create_plink:
        hl.export_plink(dense_mt, output_path, ind_id=dense_mt.s)

    b = get_batch()

    # Read in PLINK files created by Hail
    bfile = b.read_input_group(bed=f'{output_path}.bed', bim=f'{output_path}.bim', fam=f'{output_path}.fam')

    # Create GRM using GCTA
    create_GRM_j = b.new_job('Create GRM')
    create_GRM_j.image(image_path('gcta'))  # GCTA needs chromosome names to be only numbers in .bim file
    edit_chr_cmd = (
        f"awk '{{sub(/^chr/, \"\", $1); print}}' {bfile.bim} > {bfile.bim}.tmp && mv {bfile.bim}.tmp {bfile.bim}"
    )
    create_GRM_j.command(edit_chr_cmd)
    create_GRM_j.declare_resource_group(
        ofile={
            'bin': '{root}.grm.bin',
            'id': '{root}.grm.id',
            'nbin': '{root}.grm.N.bin',
        },
    )
    create_GRM_j.command(
        f'gcta --bfile {bfile} --make-grm --out {create_GRM_j.ofile}',
    )

    # Run PCA using GCTA
    run_PCA_j = b.new_job('Run PCA')
    run_PCA_j.image(image_path('gcta'))
    run_PCA_j.declare_resource_group(
        ofile={
            'eigenvec': '{root}.eigenvec',
            'eigenval': '{root}.eigenval',
        },
    )
    if relateds_to_drop := hl.read_table(relateds_to_drop):
        sgids_to_keep = set(dense_mt.s.collect()) - set(relateds_to_drop.s.collect())
        keep_file = f'{version}.indi.list'
        keep_flag = f'--keep ${{BATCH_TMPDIR}}/{keep_file}'
        # write each sg id on new line to a file
        keep_contents = '\n'.join(sgids_to_keep)
        collate_relateds_cmd = (
            f'printf "{keep_contents}" >> ${{BATCH_TMPDIR}}/{keep_file} && cat ${{BATCH_TMPDIR}}/{keep_file}'
        )
        run_PCA_j.command(collate_relateds_cmd)

    run_PCA_j.command(
        f'gcta --grm {create_GRM_j.ofile} {keep_flag if relateds_to_drop else ""} --pca {n_pcs} --out {run_PCA_j.ofile}',
    )

    b.write_output(create_GRM_j.ofile, f'{output_path}')
    b.write_output(run_PCA_j.ofile, f'{output_path}')
    b.run()


if __name__ == '__main__':
    cli_main()
