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
from cpg_utils.hail_batch import init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import get_logger
from gnomad.utils.sparse_mt import default_compute_info
from gnomad.utils.vcf import adjust_vcf_incompatible_types


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
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # not sure if we need this either
    # naive coalesce just lumps adjacent partitions together, so this could create wildly unbalanced partitions
    if partitions:
        mt = mt.repartition(partitions)

    mt.write(dense_mt_out, overwrite=True)

    if sites_only:
        mt = hl.read_matrix_table(dense_mt_out)
        info_ht = default_compute_info(mt, site_annotations=True, n_partitions=mt.n_partitions())
        info_ht = info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))
        info_ht = adjust_vcf_incompatible_types(
            info_ht,
            # with default INFO_VCF_AS_PIPE_DELIMITED_FIELDS, AS_VarDP will be converted
            # into a pipe-delimited value e.g.: VarDP=|132.1|140.2
            # which breaks VQSR parser (it doesn't recognise the delimiter and treats
            # it as an array with a single string value "|132.1|140.2", leading to
            # an IndexOutOfBound exception when trying to access value for second allele)
            pipe_delimited_annotations=[],
        )
        hl.export_vcf(info_ht, sites_only, tabix=True, parallel='header_per_shard')


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
        sites_only=args.site_only,
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()
