#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

import logging

import click

from cpg_pipes.utils import exists
from cpg_pipes.pipeline import create_pipeline, pipeline_click_options
from cpg_pipes.stages.somalier import CramSomalierAncestry, CramSomalierPedigree

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
        name='pedigree',
        stages=[
            CramSomalierPedigree,
            CramSomalierAncestry,
        ],
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
