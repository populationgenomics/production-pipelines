import logging
import struct
from argparse import ArgumentParser
from ast import Dict

import numpy as np

import hail as hl
import hailtop.batch as hb
from hailtop.batch.resource import ResourceGroup

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch


def read_grm_bin(prefix, all_n=False, size=4):
    def sum_i(i):
        return sum(range(1, i + 1))

    # File paths
    bin_file_name = prefix['grm.bin']
    n_file_name = prefix['grm.N.bin']
    id_file_name = prefix['grm.id']

    # Read ID file
    id_data = np.loadtxt(id_file_name, dtype=str)
    n = id_data.shape[0]

    # Read GRM bin file
    with open(bin_file_name, "rb") as bin_file:
        grm = np.fromfile(bin_file, dtype=np.float32 if size == 4 else np.float64, count=n * (n + 1) // 2)

    # Read N bin file
    with open(n_file_name, "rb") as n_file:
        if all_n:
            N = np.fromfile(n_file, dtype=np.float32 if size == 4 else np.float64, count=n * (n + 1) // 2)
        else:
            N = struct.unpack('f' if size == 4 else 'd', n_file.read(size))[0]

    # Compute indices for diagonal elements
    i = np.array([sum_i(j) - 1 for j in range(1, n + 1)])

    return {
        "diag": grm[i],  # Diagonal elements of GRM
        "off": np.delete(grm, i),  # Off-diagonal elements
        "id": id_data,  # ID data
        "N": N,  # N values
    }


def cli_main():

    parser = ArgumentParser()
    parser.add_argument('--output-path', help='Path to folder containingGRM files')
    parser.add_argument('--version', help='Version of the plink files')
    parser.add_argument('--n-pcs', help='Number of PCs to compute', default=10)
    parser.add_argument('--relateds-to-drop', help='HT of relateds to drop', default=None)
    args = parser.parse_args()

    run_PCA(
        output_path=args.prefix,
        version=args.version,
        n_pcs=args.n_pcs,
        relateds_to_drop=args.relateds_to_drop,
    )


def run_PCA(
    b: hb.Batch,
    grm_directory: str,
    output_path: str,
    version: str,
    n_pcs: int,
    relateds_to_drop: str,
):

    init_batch()

    # Create PCA job
    run_PCA_j = b.new_job('Run PCA')
    run_PCA_j.image(image_path('gcta'))
    # Read GRM files
    bfile = b.read_input_group(
        **{
            'grm.bin': f'{grm_directory}.grm.bin',
            'grm.id': f'{grm_directory}.grm.id',
            'grm.N.bin': f'{grm_directory}.grm.N.bin',
        },
    )
    run_PCA_j.declare_resource_group(
        ofile={
            'eigenvec': '{root}.eigenvec',
            'eigenval': '{root}.eigenval',
        },
    )

    if relateds_to_drop := hl.read_table(relateds_to_drop):
        sgids_to_remove = set(relateds_to_drop.s.collect())
        remove_file = f'{version}.indi.list'
        remove_flag = f'--remove ${{BATCH_TMPDIR}}/{remove_file}'
        # Write each sg id on new line to a file
        grm_data = read_grm_bin(bfile, all_n=True)
        remove_contents = ''
        for fam_id, sg_id in grm_data['id']:
            if sg_id in sgids_to_remove:
                remove_contents += f'{fam_id}\t{sg_id}\n'
        collate_relateds_cmd = (
            f'printf "{remove_contents}" >> ${{BATCH_TMPDIR}}/{remove_file} && cat ${{BATCH_TMPDIR}}/{remove_file}'
        )
        run_PCA_j.command(collate_relateds_cmd)

    run_PCA_j.command(
        f'gcta --grm {bfile} {remove_flag if relateds_to_drop else ""} --pca {n_pcs} --out {run_PCA_j.ofile}',
    )

    b.write_output(run_PCA_j.ofile, f'{output_path}')

    return run_PCA_j


if __name__ == '__main__':
    cli_main()
