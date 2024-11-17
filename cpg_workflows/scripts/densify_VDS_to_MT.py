#!/usr/bin/env python3

"""
takes a path to a VDS, and an output path
densifies the VDS into a MatrixTable
Writes the MT to the output path
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import config_retrieve
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import get_logger
from cpg_utils.hail_batch import init_batch


def main(vds_in: str, dense_mt_out: str, partitions: int | None = None) -> hl.MatrixTable:
    """
    Load a sparse VariantDataset
    Write a dense MatrixTable with split multiallelics
    optionally coalesce the data into a predetermined number of partitions

    Args:
        vds_in (str):
        dense_mt_out (str):
        partitions (int): if specified, write data as this many partitions
    """

    init_batch()

    # if we need to manually specify a non-standard Hail QoB JAR file
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    vds = hl.vds.read_vds(vds_in)

    get_logger(__file__).info('Splitting multiallelics')
    vds = hl.vds.split_multi(vds)

    get_logger(__file__).info('Densifying data...')
    mt = hl.vds.to_dense_mt(vds)

    # not sure if we need this either
    # naive coalesce just lumps adjacent partitions together, so this could create wildly unbalanced partitions
    if partitions:
        mt = mt.naive_coalesce(partitions)

    return mt.write(dense_mt_out, overwrite=True)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the input VDS', required=True,)
    parser.add_argument('--output', help='Path to write the result', required=True,)
    parser.add_argument('--partitions', help='Number of partitions to write the MT as', type=int,)
    args = parser.parse_args()
    main(vds_in=args.input, dense_mt_out=args.output)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()