#!/usr/bin/env python3

"""
Collect all the dataset TSVs and save them in a single file
"""

import json
from argparse import ArgumentParser

from cpg_utils import to_path


def main(mapping: str, project: str, input_file: str):
    """
    Collect all the dataset TSVs and save them in a single file
    Embellish the data with the mapping and project name

    Args:
        mapping ():
        project ():
        input_file ():

    Returns:

    """

    # take one example file, find other input files using this
    file_as_path = to_path(input_file)
    parent_dir = file_as_path.parent
    input_files = parent_dir.glob('*.tsv')

    data = {}
    for file in input_files:
        with open(file) as handle:
            data.update(json.load(handle))

    with open(output_file, 'w') as handle:
        json.dump(data, handle, indent=2)

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--mapping', help='Mapping of CPG ID to External IDs')
    parser.add_argument('--project', help='Name of Seqr Project')
    parser.add_argument('--input_example', help='A single input file')
    parser.add_argument('--output_file', help='Output file to save the combined data')
    args = parser.parse_args()

    data = {}
    for file in args.input_files:
        with open(file) as handle:
            data.update(json.load(handle))

    with open(args.output_file, 'w') as handle:
        json.dump(data, handle, indent=2)