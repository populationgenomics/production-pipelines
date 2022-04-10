#!/usr/bin/env python3

"""
Batch pipeline to generate CRAM and GVCF for all samples with alignment inputs.
"""

import logging
import click

from cpg_pipes.pipeline import pipeline_click_options, create_pipeline
from cpg_pipes.stages.align import Align
from cpg_pipes.stages.gvcf import GenotypeSample

logger = logging.getLogger(__file__)


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):
    """
    Entry point, decorated by pipeline click options.
    """
    pipeline = create_pipeline(
        name='cram_gvcf',
        description='CRAM+GVCF',
        stages=[Align, GenotypeSample],
        **kwargs,
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
