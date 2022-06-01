#!/usr/bin/env python3

"""
Batch pipeline to generate CRAMs for all samples with alignment inputs.
"""

import logging
import click

from cpg_pipes.pipeline import pipeline_options, create_pipeline
from cpg_pipes.providers.refdata import RefData
from cpg_pipes.stages.align import Align

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--cram-qc/--no-cram-qc',
    'cram_qc',
    default=True,
    is_flag=True,
    help='Run CRAM QC and PED checks',
)
@click.option(
    '--realignment-shards-num',
    'realignment_shards_num',
    type=click.INT,
    default=RefData.number_of_shards_for_realignment,
    help='Number of shards to parallelise realignment',
)
@click.option(
    '--realign',
    '--realign-from-cram-version',
    'realign_from_cram_version',
    help='Realign CRAM whenever available, instead of using FASTQ. '
         'The parameter value should correspond to CRAM version '
         '(e.g. v0 in gs://cpg-fewgenomes-main/cram/v0/CPG0123.cram)'
)
@pipeline_options
def main(**kwargs):
    """
    Entry point, decorated by pipeline click options.
    """
    kwargs['name'] = kwargs.get('name', 'cram')
    pipeline = create_pipeline(
        stages=[Align],
        **kwargs,
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
