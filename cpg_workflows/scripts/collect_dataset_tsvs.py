#!/usr/bin/env python3

"""
Collect all the dataset TSVs and save them in a single file
"""

import json
from argparse import ArgumentParser



if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('mapping', help='Mapping of CPG ID to External IDs')
    parser.add_argument('project', help='Name of Seqr Project')
    parser.add_argument('input_files', nargs='+', help='List of input files to be combined')
    parser.add_argument('output_file', help='Output file to save the combined data')
    args = parser.parse_args()

    data = {}
    for file in args.input_files:
        with open(file) as handle:
            data.update(json.load(handle))

    with open(args.output_file, 'w') as handle:
        json.dump(data, handle, indent=2)