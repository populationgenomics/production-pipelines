#!/usr/bin/python
"""
Purpose: test GCP upload from Hail Batch with multipart uploads disabled

Small script for generating a hail Batch, containing a single job.
This job will localise a file and delocalise to a new path.

This is not intended to be merged into main, it can be used sufficiently in test buckets/feature branch
"""

from argparse import ArgumentParser

from cpg_utils import hail_batch, to_path


def main(input_file: str, output_file: str):
    if to_path(output_file).exists():
        raise ValueError(f"Output file {output_file} already exists")

    batch_instance = hail_batch.get_batch('test multipart upload')

    localised_file = batch_instance.read_input(input_file)

    job = batch_instance.new_bash_job('The Job')
    job.command(f'ls {localised_file}')

    batch_instance.write_output(localised_file, output_file)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    main(input_file=args.input, output_file=args.output)
