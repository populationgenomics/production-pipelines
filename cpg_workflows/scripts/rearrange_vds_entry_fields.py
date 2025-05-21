#!/usr/bin/python3

"""
This script rearranges the fields in a VDS entry to match the order of the fields in the VDS schema.
"""

import argparse
import logging

import hail as hl


def main(modify_vds: str, match_vds: str, temp_mt: str) -> None:
    """
    Rearranges the fields in a VDS entry to match the order of the fields in the VDS schema.
    This is done by reading the reference data MatrixTable in each, reformatting to match the `match_vds` schema,
    then overwriting the `modify_vds` reference_data MatrixTable with the new schema.

    Args:
        modify_vds (str): The VDS to modify, matching it to the schema of `match_vds`
        match_vds (str): The VDS to match the schema against
        temp_mt (str): Where to write the modified MatrixTable to (can't overwrite where the MT is being read from)

    Returns:
        None
    """

    # set as a variable, we'll use this a couple of times
    modify_ref_data = f'{modify_vds}/reference_data'

    logging.info(f'Reading reference data from {modify_ref_data}')

    modify_ref_mt = hl.read_matrix_table(modify_ref_data)
    match_ref_mt = hl.read_matrix_table(f'{match_vds}/reference_data')

    logging.info('schema before modification:')
    modify_ref_mt.describe()

    # modify the reference data to match the schema of the match VDS
    modify_ref_mt = modify_ref_mt.select_entries(*match_ref_mt.entry)

    logging.info('schema after modification:')
    modify_ref_mt.describe()

    logging.info(f'Starting write to {temp_mt}')
    # write the modified reference data to a temp file
    modify_ref_mt = modify_ref_mt.checkpoint(temp_mt, overwrite=True)

    logging.info('Finished writing to temp file, overwriting initial reference data')
    # now overwrite the original reference data with the modified one
    modify_ref_mt.write(modify_ref_data, overwrite=True)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('--modify', required=True, help='Input VDS file')
    parser.add_argument('--match', required=True, help='VDS file to match schema against')
    parser.add_argument('--temp', required=True, help='A temp file to write the rearranged VDS to')
    args = parser.parse_args()

    main(modify_vds=args.modify, match_vds=args.match, temp_mt=args.temp)
