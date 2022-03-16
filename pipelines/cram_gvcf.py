#!/usr/bin/env python3

"""
Batch pipeline to generate CRAM and GVCF for all samples with alignment inputs.
"""

import logging

import click

from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.cram import CramStage
from cpg_pipes.stages.gvcf import GvcfStage

logger = logging.getLogger(__file__)


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):  # pylint: disable=missing-function-docstring

    pipeline = Pipeline(
        name='cram_gvcf',
        description='CRAM+GVCF',
        stages_in_order=[CramStage, GvcfStage],
        **kwargs,
    )
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
