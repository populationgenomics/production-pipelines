#! /usr/bin/env python

"""
This script is used to convert a hail MatrixTable/Table/VDS to a vcf file.
Prototype design using hail's parallel export_vcf method.

Doubles up by using the gcloud compose command to concatenate the VCFs
without ever needing to localise the separate elements.

Usage:
    make_vcf_from_hail.py <input> <output> [--overwrite] [--sites_only] [--temp <gs://temp>]

Plan:

1. Load the hail object
2. Export to vcf using hail's export_vcf method in parallel
3. parse the manifest file to identify the chunks in order
4. Generate a series of Bash jobs to cat together chunks of the output (remove?)
  - use the GCP entity stat to estimate storage requirements
5. Concatenate the vcf intermediates
6. Index the concatenated vcf file using tabix
7. Write combined file to GCP
"""

import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, image_path, init_batch
from cpg_workflows.utils import chunks


LOGGER = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
LOGGER.addHandler(handler)
LOGGER.setLevel(logging.INFO)

# make this configurable
MIN_NUM_PARTITIONS = 100


def read_hail(path: str) -> hl.Table | hl.MatrixTable | hl.vds.VariantDataset:
    """
    read a hail object using the appropriate method
    if a MT or Table, require a minimum number of partitions to
    make this parallel operation worthwhile
    TODO is this repartitioning valid?

    Args:
        path (str): path to the input object
    Returns:
        hail object (hl.VDS, hl.MatrixTable, or hl.Table), and a data type
    """
    LOGGER.info(f'Reading data from {path}')
    if path.strip('/').endswith('.vds'):
        # not sure if a VDS has a count/partitions equivalent
        LOGGER.info('Reading VDS')
        t = hl.vds.read_vds(path)
        return t

    if path.strip('/').endswith('.ht'):
        LOGGER.info('Reading Table')
        t = hl.read_table(path)
        t = hl.read_matrix_table(
            path, _n_partitions=max(t.n_partitions(), MIN_NUM_PARTITIONS)
        )
    else:
        assert path.strip('/').endswith('.mt')
        LOGGER.info('Reading MT')
        t = hl.read_matrix_table(path)
        t = hl.read_matrix_table(
            path, _n_partitions=max(t.n_partitions(), MIN_NUM_PARTITIONS)
        )
    LOGGER.info(
        f'Read data from {path} - dimensions: {t.count()}, partitions: {t.n_partitions()}'
    )
    return t


def vds_processing(
    obj: hl.vds.VariantDataset, sites_only: bool = False
) -> hl.MatrixTable | hl.Table:
    """
    Process a VDS object - either densify if we need the full MT
    or just break off the variant MT.rows() if sites_only

    # needs more stuff here, there's some incompatible VDS fields

    Args:
        obj (hl.vds.VariantDataset): hail object
        sites_only (bool): Remove genotypes from the VCF representation

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


def get_gcloud_condense_command_string(fragment_list: list[str], output: str) -> str:
    """
    Generate a gcloud compose command string to concatenate a list of VCF fragments
    write the product to the output path

    Args:
        fragment_list ():
        output ():

    Returns:

    """

    cat_string = 'gcloud storage objects compose '
    for fragment in fragment_list:
        cat_string += f'{fragment} '
    cat_string += f'{output}'
    return cat_string


def compose_condense_fragments(
        path_list: list[str],
        temp_prefix: str | None = None,
        chunk_size: int = 30,
        output_name: str | None = None
) -> list[str]:
    """
    takes a list of things to condense, creates jobs to condense them
    and returns a list of the result paths

    Args:
        path_list ():
        temp_prefix ():
        chunk_size (): number of objects to condense at once
        output_name (str): if specified, write

    Returns:
        list of paths to the condensed objects
    """

    assert output_name or temp_prefix, 'Either output_name or temp_prefix must be specified'

    chunk_job = get_batch().new_job(name=f'chunkbuster_{len(path_list)}')
    chunk_job.image(get_config()['workflow']['driver_image'])
    authenticate_cloud_credentials_in_job(chunk_job)

    sub_chunks = []
    for i, chunk in enumerate(chunks(path_list, chunk_size=chunk_size)):
        # either the named file, or a temp location
        chunk_output = output_name or join(temp_prefix, f'chunk_{i}.vcf.bgz')
        chunk_job.command(get_gcloud_condense_command_string(chunk, chunk_output))
        sub_chunks.append(chunk_output)
    return sub_chunks


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
    if temp is None:
        temp = get_config()['storage']['default']['tmp']
        LOGGER.info(f'Using default temp path {temp}')

    assert isinstance(temp, str), 'Temp path must be a string'
    assert temp.startswith('gs://'), 'Temp must be a GCP path'
    LOGGER.info(f'Using temp path {temp}')

    # temp needs a `.vcf.bgz` extension, otherwise data is dumped in plain text
    if not temp.endswith('vcf.bgz'):
        temp = join(temp, 'temp.vcf.bgz')
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

    # fourth - parse the manifest file and generate commands to concatenate vcf files
    # done in-line here, but equally this could be a separate (python?) job
    # read lines, create a script
    # next job localises the script and runs
    # this current approach wouldn't work within a pipeline
    LOGGER.info(
        'Parsing the manifest file and generating a bash script for concatenating the VCF files'
    )
    manifest = hl.hadoop_open(join(temp, 'shard-manifest.txt')).readlines()

    # prefix these shard names to get full GCP path for each
    chunk_paths = [join(temp, str(fragment).strip()) for fragment in manifest]

    # rolling squash of the chunks, should enable scaling past 900 shards (2 rounds)
    temp_chunk_prefix_num = 1
    while len(chunk_paths) > 30:
        condense_temp = join(temp, f'temp_chunk_{temp_chunk_prefix_num}')
        temp_chunk_prefix_num += 1
        chunk_paths = compose_condense_fragments(chunk_paths, condense_temp)

    # now compress all those chunks into the final output file
    final_path = compose_condense_fragments(chunk_paths, temp, output_name=output_path)[0]

    # todo - after completion remove the temp files/dir?

    # one final job - read the final vcf in, index it, move the index, and write it out
    input_vcf = get_batch().read_input(final_path)
    final_job = get_batch().new_bash_job(name='index_final_vcf')
    final_job.image(image_path('bcftools'))
    final_job.command(f'tabix {input_vcf} && mv {input_vcf}.tbi {final_job.index}')
    get_batch().write_output(final_job.index, final_path + '.tbi')
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
