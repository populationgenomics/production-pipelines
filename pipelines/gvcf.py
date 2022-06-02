#!/usr/bin/env python3

"""
Batch pipeline to generate GVCF for all samples CRAMs.
"""

import logging
import click

from cpg_pipes.jobs import haplotype_caller
from cpg_pipes.pipeline import pipeline_options, create_pipeline
from cpg_pipes.stages.genotype_sample import GenotypeSample
from cpg_pipes.utils import exists

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--hc-intervals-num',
    'hc_intervals_num',
    type=click.INT,
    default=haplotype_caller.DEFAULT_INTERVALS_NUM,
    help='Number of intervals to devide the genome for sample genotyping with '
    'gatk HaplotypeCaller',
)
@pipeline_options
def main(
    **kwargs,
):
    """
    Entry point, decorated by pipeline click options.
    """
    kwargs['name'] = kwargs.get('name', 'gvcf')
    pipeline = create_pipeline(
        stages=[GenotypeSample],
        **kwargs,
    )
    if pipeline.skip_samples_with_missing_input:
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(f'Could not find CRAM, skipping sample {sample.id}')
                sample.active = False

    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
