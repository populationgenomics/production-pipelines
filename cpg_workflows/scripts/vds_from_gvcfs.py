#! /usr/bin/env python3

"""
Takes a bunch of gVCFs, makes a VDS
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, output_path
from cpg_workflows.utils import get_logger


def extract_vcf_from_dataset_vcf(gvcfs: list[str], sgids: list[str], out_path: str, tmp_prefix: str | None = None):
    """
    Takes paths to a bunch of gVCFs
    Combines into a VDS
    writes that out
    Args:
        gvcfs ():
        sgids ():
        out_path ():
        tmp_prefix ():
    """

    if tmp_prefix is None:
        tmp_prefix = output_path('vds_combiner', category='tmp')
    if get_config()['workflow']['sequencing_type'] == 'exome':
        params = {'use_exome_default_intervals': True}
    else:
        params = {'use_genome_default_intervals': True}

    combiner = hl.vds.new_combiner(
        gvcf_paths=gvcfs,
        gvcf_sample_names=sgids,
        # Header must be used with gvcf_sample_names, otherwise gvcf_sample_names
        # will be ignored. The first gvcf path works fine as a header because it will
        # be only read until the last line that begins with "#":
        gvcf_external_header=gvcfs[0],
        output_path=str(out_path),
        reference_genome='GRCh38',
        temp_path=tmp_prefix,
        force=True,
        **params,
    )
    combiner.run()


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--gvcfs', help='Input gVCFs', nargs='+')
    parser.add_argument('--sgids', help='Input SG IDs', nargs='+')
    parser.add_argument('--out', help='Output VCF path')
    args = parser.parse_args()
    assert len(args.gvcfs) == len(args.sgids)
    get_logger(__file__).info(f'Creating VDS {args.out} from {len(args.gvcfs)} gVCFs')

    init_batch()
    extract_vcf_from_dataset_vcf(gvcfs=args.gvcfs, sgids=args.sgids, out_path=args.out)
