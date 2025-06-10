"""
Hail 0.2.134 is highly problematic - it has the glaring limitation that it cannot combine data into any VDS created
before Hail 134, which at the moment is almost all our data.

This script takes a path to a VDS on 133, writes a copy to a given location. That location can be passed to a new
combiner run, which will then be able to combine data into it.
"""

import argparse

import hail as hl

from cpg_utils import hail_batch

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Write Hail 0.2.133 VDS as 0.2.134, so that it can be used to combine')
    parser.add_argument('vds_path', type=str, help='Path to the VDS to swap out')
    parser.add_argument('tmp_path', type=str, help='Temporary path to write the swapped VDS to')
    args = parser.parse_args()

    hail_batch.init_batch()

    vds = hl.vds.read_vds(args.vds_path)
    vds.write(args.tmp_path, overwrite=True)
