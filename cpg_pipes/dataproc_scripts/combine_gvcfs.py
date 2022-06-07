#!/usr/bin/env python

"""
Combine a set of GVCFs into a MatrixTable
"""

import os
import subprocess
from os.path import join, basename
import logging
import shutil
from pathlib import Path

import pandas as pd
import click
import hail as hl
from cpg_utils.hail_batch import genome_build

from cpg_pipes.hailquery import init_hail
from cpg_pipes import version
from cpg_pipes.pipeline.entry import file_validation_callback

logger = logging.getLogger(__name__)


@click.command()
@click.version_option(version.__version__)
@click.option(
    '--meta-csv',
    'meta_csv_path',
    required=True,
    callback=file_validation_callback(ext='csv', must_exist=True),
    help='Sample data CSV path',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    required=True,
    callback=file_validation_callback(ext='mt'),
    help='path to write the combined MatrixTable',
)
@click.option(
    '--bucket',
    'work_bucket',
    required=True,
    help='path to folder for intermediate output. '
    'Can be a Google Storage URL (i.e. start with `gs://`).',
)
@click.option(
    '--local-tmp-dir',
    'local_tmp_dir',
    help='local directory for temporary files and Hail logs (must be local).',
)
def main(
    meta_csv_path: str,
    out_mt_path: str,
    work_bucket: str,
    local_tmp_dir: str,
):  # pylint: disable=missing-function-docstring
    """
    Entry point
    """
    init_hail('combine_gvcfs', Path(local_tmp_dir))

    # Copy the metadata file locally
    local_meta_csv_path = join(local_tmp_dir, basename(meta_csv_path))
    gsutil_cp(meta_csv_path, local_meta_csv_path)
    samples_df = pd.read_table(local_meta_csv_path)

    logger.info(f'Combining {len(samples_df.gvcf_path)} GVCFs')
    hl.experimental.run_combiner(
        list(samples_df.gvcf_path),
        sample_names=list(samples_df.s),  # pylint: disable=no-member
        out_file=out_mt_path,
        reference_genome=genome_build(),
        use_genome_default_intervals=True,
        tmp_path=os.path.join(work_bucket, 'tmp'),
        overwrite=True,
        key_by_locus_and_alleles=True,
    )

    mt = hl.read_matrix_table(out_mt_path)
    logger.info(
        f'Written {mt.cols().count()} samples to {out_mt_path}, '
        f'n_partitions={mt.n_partitions()}'
    )

    shutil.rmtree(local_tmp_dir)


def gsutil_cp(
    src_path: str,
    dst_path: str,
    disable_check_hashes: bool = False,
    recursive: bool = False,
    quiet: bool = False,
) -> str:
    """
    Wrapper around `gsutil cp`

    @param src_path: path to a file to copy from
    @param dst_path: path to copy to
    @param disable_check_hashes:
        Uses the gsutil option `-o GSUtil:check_hashes=never` which is required to
        get around the gsutil integrity checking error, as conda gsutil doesn't use
        CRC32c:
        > Downloading this composite object requires integrity checking with CRC32c,
          but your crcmod installation isn't using the module's C extension, so the
          hash computation will likely throttle download performance.

          To download regardless of crcmod performance or to skip slow integrity
          checks, see the "check_hashes" option in your boto config file.
    @param recursive: to copy a directory
    @param quiet: disable logging of commands and copied files
    @returns: dst_path
    """
    cmd = (
        'gsutil '
        + ('-q ' if quiet else '')
        + ('-o GSUtil:check_hashes=never ' if disable_check_hashes else '')
        + 'cp '
        + ('-r ' if recursive else '')
        + f'{src_path} {dst_path}'
    )
    if not quiet:
        logger.info(cmd)
    subprocess.run(cmd, check=False, shell=True)
    return dst_path


if __name__ == '__main__':
    main()  # pylint: disable=E1120
