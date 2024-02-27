#! /usr/bin/env python

"""
This script is used to convert a hail MatrixTable/Table/VDS to a vcf file.
Comparison design doing nothing special, just a simple conversion.
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch


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

    Args:
        obj (hl.vds.VariantDataset): hail object
        sites_only (bool): remove genotypes if appropriate

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


def create_vcf_from_hail(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    sites_only: bool = False
):
    """
    This function is used to convert a hail MatrixTable/Table/VDS to a vcf file.

    Args:
        input_path (str): Input hail object
        output_path (str): Output vcf file
        overwrite (bool): Overwrite output file if it exists
        sites_only (bool): Remove genotypes if appropriate
    """

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
    hl.export_vcf(obj, output=output_path)


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
    args = parser.parse_args()

    create_vcf_from_hail(
        args.input, args.output, args.overwrite, args.sites_only
    )
