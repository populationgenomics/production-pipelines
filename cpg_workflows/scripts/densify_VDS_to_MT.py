#!/usr/bin/env python3

"""
takes a path to a VDS, and an output path
densifies the VDS into a MatrixTable
Writes the MT to the output path
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import get_logger


def main(
    vds_in: str,
    dense_mt_out: str,
    partitions: int | None = None,
    sites_only: str | None = None,
) -> None:
    """
    Load a sparse VariantDataset
    Write a dense MatrixTable with split multiallelics
    optionally coalesce the data into a predetermined number of partitions

    Args:
        vds_in (str):
        dense_mt_out (str):
        partitions (int): if specified, write data as this many partitions
        sites_only (str): optional, if used write a sites-only VCF directory to this location
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

    # remove any monoallelic or non-ref-in-any-sample sites
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.GT.is_non_ref())))

    # naive coalesce just lumps adjacent partitions together, so this could create wildly unbalanced partitions
    # if this is not used we'll retain the original partitions
    if partitions:
        mt = mt.naive_coalesce(partitions)

    # by default, drop the GVCF info
    if "gvcf_info" in mt.entry:
        mt = mt.drop('gvcf_info')

    mt.write(dense_mt_out, overwrite=True)

    if sites_only:
        mt = hl.read_matrix_table(dense_mt_out)
        hl.export_vcf(mt, sites_only, tabix=True, parallel='header_per_shard')


def cli_main():
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to the input VDS',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Path to write the result',
        required=True,
    )
    parser.add_argument(
        '--partitions',
        help='Number of partitions to write the MT as',
        type=int,
    )
    parser.add_argument(
        '--sites_only',
        help='Specify an output path for a sites-only VCF, or None',
    )
    args = parser.parse_args()
    main(
        vds_in=args.input,
        dense_mt_out=args.output,
        partitions=args.partitions,
        sites_only=args.sites_only,
    )


if __name__ == '__main__':
    cli_main()
