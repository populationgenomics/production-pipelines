#!/usr/bin/python3

"""
This script rearranges the fields in a VDS entry to match the order of the fields in the VDS schema.

example invocation, resolving our current problem VDSs:
analysis-runner \
    --dataset seqr \
    --description 'overwrite reference VDS schema' \
    -o schema_correction \
    --access-level full \
    python3 cpg_workflows/scripts/rearrange_vds_entry_fields.py \
        --modify gs://cpg-seqr-main-tmp/rd_combiner/fa55fe991ca12256081ef2a9cf3a1218720db9_5734/CreateVdsFromGvcfsWithHailCombiner/temp_dir/combiner_removal_temp.vds \
        --match gs://cpg-seqr-main-tmp/rd_combiner/fa55fe991ca12256081ef2a9cf3a1218720db9_5734/CreateVdsFromGvcfsWithHailCombiner/temp_dir/combiner-intermediates/cdbe20b6-a1b
c-45f4-b80a-ff44578f8236_gvcf-combine_job1/dataset_0.vds \
        --temp gs://cpg-seqr-main-tmp/oneshot_patch/fa55fe991ca12256081ef2a9cf3a1218720db9_5734.mt

"""

import argparse
import logging

import hail as hl

from cpg_utils.hail_batch import init_batch


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

    init_batch()

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
