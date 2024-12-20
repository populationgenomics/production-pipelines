"""
See https://batch.hail.populationgenomics.org.au/batches/556897/jobs/1

This script is used to drop fields from the entries of a VDS, and then
overwrite the original VDS with the new entries.

This aims to solve an issue with Hail gVCF Combining, where we do hierarchical merge of

gVCFs -> VDSs

then

VDSs -> single VDS

The gVCFs -> VDS stage takes a union of all entry fields across all component gVCFs, so if a single
gVCF has an unusual field, the VDS will have that field in the schema. Combining VDS intermediates
requires all input VDSs to have an identical schema, so we have to groom out the dodgy entry fields.
Due to the way data is combined, the entries will always be in the `entry.gvcf_info` field.

Due to hail's read-from-file processing, we can't patch in place, we have to write out a new VDS, then
overwrite the original VDS with the new one.
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to input VDS',
        required=True,
    )
    parser.add_argument('--temp_location', help='Path to temporary location', required=True)
    parser.add_argument(
        '--fields',
        help='List of strings, fields to remove. So far only DRAGstr issues have been observed',
        nargs='+',
        default=['DRAGstrInfo', 'DRAGstrParams'],
    )
    args = parser.parse_args()

    logging.info('Starting VDS entry field drop and overwrite')
    logging.info(f'Input VDS: {args.input}')
    logging.info(f'Temp location: {args.temp_location}')
    logging.info(f'Fields to drop: {args.fields}')

    if not args.fields:
        raise ValueError('No fields to drop')

    init_batch()

    # read in the input VDS
    vds = hl.vds.read_vds(args.input)
    vds.variant_data.describe()
    vds.variant_data = vds.variant_data.annotate_entries(
        gvcf_info=vds.variant_data.gvcf_info.drop(
            *args.fields,
        ),
    )
    logging.info(f'Dropped fields: {args.fields}')
    vds.variant_data.describe()

    logging.info(f'Checkpointing VDS to {args.temp_location}')
    vds = vds.checkpoint(args.temp_location)

    logging.info(f'Checkpointed VDS to {args.temp_location}, writing back to {args.input}')
    vds.write(args.input, overwrite=True)

    logging.info('Done')
