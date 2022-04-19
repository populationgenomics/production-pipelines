#!/usr/bin/env python3

"""
Batch pipeline to run QC. Output is a per-sample MultiQC report.
"""

import logging

import click

from cpg_pipes.pipeline import (
    pipeline_click_options,
    create_pipeline,
)
from cpg_pipes.stages.multiqc import MultiQC
from cpg_pipes.utils import exists

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
        name='qc',
        stages=[MultiQC],
        **kwargs,
    )
    if pipeline.skip_samples_with_missing_input:
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(f'Could not find CRAM, skipping sample {sample.id}')
                sample.active = False

    pipeline.run(force_all_implicit_stages=True)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
