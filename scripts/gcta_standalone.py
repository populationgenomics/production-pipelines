import logging
import struct

import click
import numpy as np

import hail as hl
import hailtop.batch as hb

from cpg_utils import Path, to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import Batch, get_batch, init_batch
from cpg_workflows.utils import can_reuse


def make_plink(mt: hl.MatrixTable, plink_output_path: Path):
    plink_paths = [f'{plink_output_path}.{ext}' for ext in ['bed', 'bim', 'fam']]
    logging.info(f'plink_paths: {plink_paths}')
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

    grm_paths = [f'{output_path}.{ext}' for ext in ['grm.bin', 'grm.id', 'grm.N.bin']]
    logging.info(f'GRM paths: {grm_paths}')
    if can_reuse(grm_paths):
        logging.info('GRM files already exist')
        return True

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
    logging.info(f'create_GRM_j.ofile: {create_GRM_j.ofile}')
    create_GRM_j.command(
        f'gcta --bfile {bfile} --make-grm --out {create_GRM_j.ofile}',
    )

    return create_GRM_j


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


def run_PCA(
    b: hb.Batch,
    create_GRM_j: hb.ResourceGroup,
    grm_directory: str,
    version: str,
    n_pcs: int,
    relateds_to_drop: str,
):

    logging.info(
        f'create_GRM_j.ofile: {create_GRM_j.ofile["grm.bin"], create_GRM_j.ofile["grm.id"], create_GRM_j.ofile["grm.N.bin"]}',
    )

    # Create PCA job
    run_PCA_j = b.new_job('Run PCA')
    run_PCA_j.image(image_path('gcta'))
    # Read GRM files
    # Turn grm_directory into a dictionary with correct keys and values
    run_PCA_j.declare_resource_group(
        grm_directory={
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
    logging.info(
        f'{run_PCA_j.grm_directory}'
        f'{run_PCA_j.grm_directory["grm.bin"]}'
        f'{run_PCA_j.grm_directory["grm.id"]}'
        f'{run_PCA_j.grm_directory["grm.N.bin"]}',
    )
    # Check if there are relateds to drop
    if relateds_to_drop:
        with to_path(relateds_to_drop).open('r') as f:
            sgids_to_remove = set(line.strip() for line in f)
        remove_file = f'{version}.indi.list'
        remove_flag = f'--remove ${{BATCH_TMPDIR}}/{remove_file}'
        id_data = np.loadtxt(run_PCA_j.grm_directory['grm.id'], dtype=str)
        remove_contents = ''
        for fam_id, sg_id in id_data:
            if sg_id in sgids_to_remove:
                remove_contents += f'{fam_id}\t{sg_id}\n'
        collate_relateds_cmd = (
            f'printf "{remove_contents}" >> ${{BATCH_TMPDIR}}/{remove_file} && cat ${{BATCH_TMPDIR}}/{remove_file}'
        )
        run_PCA_j.command(collate_relateds_cmd)

    run_PCA_j.command(
        f'gcta --grm {run_PCA_j.grm_directory} {remove_flag if relateds_to_drop else ""} --pca {n_pcs} --out {run_PCA_j.ofile}',
    )

    return run_PCA_j


@click.command()
@click.option('--dense-mt-path')
@click.option('--out-dir')
@click.option('--relateds-to-drop-path')
@click.option('--version')
@click.option('--create-plink', is_flag=True, help="Create PLINK output based on out_dir")
@click.option('--create-grm', is_flag=True, help="Create GRM output based on out_dir")
def main(
    dense_mt_path: Path,
    out_dir: Path,
    relateds_to_drop_path: str,
    version: str,
    create_plink: bool,
    create_grm: bool,
):
    logging.basicConfig(level=logging.INFO)
    out_dir = to_path(out_dir)
    plink_output_path = out_dir / 'plink' / version
    grm_output_path: Path = out_dir / 'grm' / version

    b = get_batch()
    init_batch()
    # Make Plink files
    logging.info(f'plink_output_path: {plink_output_path}')

    logging.info(f'Loading dense MT from {dense_mt_path}')
    dense_mt = hl.read_matrix_table(dense_mt_path)
    logging.info('Loaded dense MT')
    plink_result = make_plink(dense_mt, to_path(plink_output_path))
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
            output_path=grm_output_path,
        )
        b.write_output(create_GRM_j.ofile, str(grm_output_path))

    # Run PCA
    pca_output_path = out_dir / 'pca' / version
    if relateds_to_drop_path:
        logging.info('Reading relateds_to_drop Hail Table')
        relateds_to_drop_ht = hl.read_table(relateds_to_drop_path)
        logging.info('Read relateds_to_drop, now collecting list of sample ids')
        sample_ids = relateds_to_drop_ht.s.collect()
        logging.info(f'sample ids to drop: {sample_ids}')
        relateds_txt_path = out_dir / version / 'gcta_relateds_remove.indi.list'
        logging.info(f'relateds_txt_path: {relateds_txt_path}')
        with to_path(relateds_txt_path).open('w') as f:
            for sample_id in sample_ids:
                f.write(f'{sample_id}\n')
        logging.info('Creating PCA job')
        run_PCA_j = run_PCA(
            b=b,
            create_GRM_j=create_GRM_j,
            grm_directory=grm_output_path,
            version=version,
            n_pcs=10,
            relateds_to_drop=relateds_txt_path,
        )
        run_PCA_j.depends_on(create_GRM_j)
        b.write_output(run_PCA_j.ofile, str(pca_output_path))

    # b = get_batch()
    logging.info('Running batch')
    b.run()


if __name__ == '__main__':
    main()
