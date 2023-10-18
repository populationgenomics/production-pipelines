#!/usr/bin/env python3

"""
Query script to parse JSON VEP results.
"""

import click
import hail as hl

from cpg_utils.hail_batch import genome_build
from cpg_workflows.query_modules import vep


@click.command()
@click.argument('vep_results_paths', nargs=-1)
@click.option('--out_path', required=True)
@click.option('--vep_version', required=True)
def main(vep_results_paths: list[str], out_path: str, vep_version: str):
    hl.init(default_reference=genome_build())

    vep.vep_json_to_ht(
        vep_result_paths=vep_results_paths,
        out_path=out_path,
        vep_version=vep_version
    )


if __name__ == '__main__':
    main()
