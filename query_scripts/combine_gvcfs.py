#!/usr/bin/env python

"""
Combine a set of GVCFs into a MatrixTable
"""

import os
from os.path import join, basename
import logging
import shutil
from pathlib import Path

import pandas as pd
import click
import hail as hl

from cpg_pipes.hailquery import init_hail
from cpg_pipes.utils import get_validation_callback
from cpg_pipes import buckets, utils, _version

logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


@click.command()
@click.version_option(_version.__version__)
@click.option(
    '--meta-csv',
    'meta_csv_path',
    required=True,
    callback=get_validation_callback(ext='csv', must_exist=True),
    help='Sample data CSV path',
)
@click.option(
    '--out-mt',
    'out_mt_path',
    required=True,
    callback=get_validation_callback(ext='mt'),
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
@click.option(
    '--hail-billing',
    'hail_billing',
    help='Hail billing account ID.',
)
def main(
    meta_csv_path: str,
    out_mt_path: str,
    work_bucket: str,
    local_tmp_dir: str,
    hail_billing: str,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring
    local_tmp_dir = init_hail('combine_gvcfs', Path(local_tmp_dir))

    # Copy the metadata file locally    
    local_meta_csv_path = join(local_tmp_dir, basename(meta_csv_path))
    buckets.gsutil_cp(meta_csv_path, local_meta_csv_path)
    samples_df = pd.read_table(local_meta_csv_path)

    logger.info(f'Combining {len(samples_df.gvcf_path)} GVCFs')
    hl.experimental.run_combiner(
        list(samples_df.gvcf_path),
        sample_names=list(samples_df.s),  # pylint: disable=no-member
        out_file=out_mt_path,
        reference_genome=utils.DEFAULT_REF,
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


if __name__ == '__main__':
    main()  # pylint: disable=E1120
