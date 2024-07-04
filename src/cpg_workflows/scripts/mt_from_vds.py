#! /usr/bin/env python3

"""
Takes a bunch of gVCFs, makes a VDS
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--vds', help='Input VDS')
    parser.add_argument('--mt', help='Output MT')
    args = parser.parse_args()

    get_logger(__file__).info(f'Creating MT {args.mt} from VDS {args.vds}')

    init_batch()

    vds = hl.vds.read_vds(args.vds)
    vds = hl.vds.split_multi(vds)
    mt = hl.vds.to_dense_mt(vds)
    mt.write(args.mt)
