#! /usr/bin/env python

"""
This script is used to convert a hail MatrixTable/Table/VDS to a vcf file.
Prototype design using hail's parallel export_vcf method.

Usage:
    make_vcf_from_hail.py <input> <output> [--overwrite]

Plan:

1. Load the hail object
2. Export to vcf using hail's export_vcf method in parallel
3. Python Job to parse the manifest file and generate a bash script for concatenating the vcf files
4. Run the bash script to concatenate the vcf files (assuming bgzip'd output files can be directly concatenated)
5. Index the concatenated vcf file using tabix
6. Write combined file to GCP
"""

import logging
import os.path
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, image_path
from cpg_workflows.utils import chunks


LOGGER = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
handler.setFormatter(formatter)
LOGGER.addHandler(handler)
LOGGER.setLevel(logging.INFO)

# make this configurable
MIN_NUM_PARTITIONS = 100


def read_hail(path: str) -> hl.Table | hl.MatrixTable | hl.vds.VariantDataset:
    """
    read a hail object using the appropriate method
    TODO is this repartitioning valid?

    Args:
        path (str): path to the input object
    Returns:
        hail object (hl.VDS, hl.MatrixTable, or hl.Table), and a data type
    """
    LOGGER.info(f'Reading data from {path}')
    if path.strip('/').endswith('.ht'):
        LOGGER.info('Reading Table')
        t = hl.read_table(path)
        t = hl.read_matrix_table(path, _n_partitions=max(t.n_partitions(), MIN_NUM_PARTITIONS))
    elif path.strip('/').endswith('.vds'):
        LOGGER.info('Reading VDS')
        t = hl.vds.read_vds(path)
    else:
        assert path.strip('/').endswith('.mt')
        LOGGER.info('Reading MT')
        t = hl.read_matrix_table(path)
        t = hl.read_matrix_table(path, _n_partitions=max(t.n_partitions(), MIN_NUM_PARTITIONS))
    LOGGER.info(f'Read data from {path} - dimensions: {t.count()}, partitions: {t.n_partitions()}')
    return t


def vds_processing(obj: hl.vds.VariantDataset, sites_only: bool = False) -> hl.MatrixTable | hl.Table:
    """
    Process a VDS object - either densify if we need the full MT
    or just break off the variant MT.rows() if sites_only

    Args:
        obj (hl.vds.VariantDataset): hail object

    Returns:
        hail object
    """
    LOGGER.info('Processing VDS into MT')
    # don't bother densifying, gimme those variants
    if sites_only:
        return obj.variant_data.rows()
    else:
        return hl.vds.to_dense_mt(obj)


def make_sites_only(obj: hl.MatrixTable | hl.Table) -> hl.MatrixTable | hl.Table:
    """
    Remove genotypes from the hail object if appropriate

    Args:
        obj (hl.MatrixTable | hl.Table | hl.vds.VariantDataset): hail object
    Returns:
        hail object with genotypes removed if appropriate
    """
    LOGGER.info('Removing genotypes if appropriate')
    if isinstance(obj, hl.MatrixTable):
        return obj.rows()
    # if isinstance(obj, hl.vds.VariantDataset):
    #     return obj.variant_data.rows()
    if isinstance(obj, hl.Table):
        return obj
    LOGGER.info('No further processing applied, assuming VDS')
    return obj

    raise ValueError(f'Unexpected hail object type {type(obj)}')


def create_vcf_from_hail(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    sites_only: bool = False,
    temp: str | None = None,
):
    """
    This function is used to convert a hail MatrixTable/Table/VDS to a vcf file.

    Args:
        input_path (str): Input hail object
        output_path (str): Output vcf file
        overwrite (bool): Overwrite output file if it exists
        sites_only (bool): Remove genotypes if appropriate
        temp (str): temp path for partial VCFs, optional
    """
    if temp:
        assert temp.startswith('gs://'), 'Temp path must be a GCP path'
        LOGGER.info(f'Using temp path {temp}')
    else:
        temp = get_config()['storage']['default']['tmp']
        LOGGER.info(f'Using default temp path {temp}')

    # temp needs a `.vcf.bgz` extension, otherwise data is dumped in plain text
    if not temp.endswith('vcf.bgz'):
        temp = os.path.join(temp, 'temp.vcf.bgz')
        LOGGER.info(f'Using edited temp path {temp}')

    # first - decide what to do?
    assert output_path.endswith('vcf.bgz'), 'Output file must end with vcf.bgz'
    if not overwrite and to_path(output_path).exists():
        raise FileExistsError(f'Output file {output_path} already exists')

    # second - load the hail object
    init_batch()
    obj = read_hail(input_path)

    # some object-specific processing
    if isinstance(obj, hl.vds.VariantDataset):
        obj = vds_processing(obj, sites_only)
    else:
        # do some stuff here if sites-only
        obj = make_sites_only(obj) if sites_only else obj

    # do any filtering/re-partitioning here

    # third - export to vcf using hail's export_vcf method in parallel
    LOGGER.info('Exporting to multiple temp VCFs')
    hl.export_vcf(obj, output=temp, parallel='separate_header')

    # fourth - parse the manifest file and generate a bash script for concatenating the vcf files
    # done in-line here, but equally this could be a separate (python?) job
    # read lines, create a script
    # next job localises the script and runs
    # this current approach wouldn't work within a pipeline
    LOGGER.info('Parsing the manifest file and generating a bash script for concatenating the VCF files')
    manifest = hl.hadoop_open(join(temp, 'shard-manifest.txt')).readlines()


    # there's a ton of possible approaches here - like doing a rolling merge
    # to reduce the overall amount of space, or splitting the merge into a bunch of different jobs.
    # this is an overly cautious approach, just to see what happens
    total_gb = 0
    sub_chunks = []
    for i, chunk in enumerate(chunks(manifest, 20)):
        LOGGER.info(f'Processing chunk {chunk}')
        chunk_job = get_batch().new_job(name=f'chunk_{i}')
        # estimate space consumed by this chunk
        storage_bytes = 0

        cat_string = 'cat '
        for fragment in chunk:
            full_frag_path = join(temp, str(fragment).strip())
            LOGGER.info(f'Processing fragment {full_frag_path}')
            storage_bytes += to_path(join(temp, full_frag_path)).stat().st_size
            frag_local = chunk_job.read_input(join(temp, full_frag_path))
            cat_string += f' {frag_local} '
        cat_string += f' > {chunk_job.output}'
        sub_chunks.append(chunk_job.output)
        chunk_job.command(cat_string)

        # total, plus margin for error, then doubled for new output file
        storage_gb = (storage_bytes * 2.2) // 1024 ** 3
        total_gb += storage_gb
        chunk_job.storage(storage_gb)

    # and now one big job
    final_job = get_batch().new_job(name='final_merge')
    final_job.image(image_path('bcftools'))
    final_job.declare_resource_group(output={'vcf.bgz': '{root}', 'vcf.bgz.tbi': '{root}.tbi'})
    command = 'cat '
    for sub_chunk in sub_chunks:
        command += f' {sub_chunk} '
    command += f' > {final_job.output}'
    final_job.storage(total_gb)
    final_job.command(command)
    final_job.command(f'tabix {final_job.output}')

    get_batch().write_output(final_job.output, output_path.removesuffix('vcf.bgz'))
    get_batch().run(wait=False)


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('input', help='Input hail object')
    parser.add_argument('output', help='Output vcf file')
    parser.add_argument(
        '--overwrite', action='store_true', help='Overwrite output file if it exists'
    )
    parser.add_argument(
        '--sites_only', action='store_true', help='Remove genotypes if appropriate'
    )
    parser.add_argument(
        '--temp', help='temp path for partial VCFs, optional', required=False
    )
    args = parser.parse_args()

    create_vcf_from_hail(
        args.input, args.output, args.overwrite, args.sites_only, args.temp
    )
