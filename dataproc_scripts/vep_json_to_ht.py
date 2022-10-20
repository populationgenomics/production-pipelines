#!/usr/bin/env python3

"""
Query script to parse JSON VEP results.
"""

import click
import hail as hl

from cpg_utils.hail_batch import genome_build
from query_modules import vep


@click.command()
@click.option(
    '--vep-result-json-path',
    'vep_results_paths',
    multiple=True,
    required=True,
)
@click.option(
    '--out-path',
    'out_path',
    required=True,
)
def main(
    vep_results_paths: list[str],
    out_path: str,
):
    hl.init(default_reference=genome_build())

    vep.vep_json_to_ht(
        vep_results_paths=vep_results_paths,
        out_path=out_path,
    )
