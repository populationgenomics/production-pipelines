#!/usr/bin/env python3

"""
Batch pipeline to generate GVCF for all samples CRAMs.
"""

import logging

from cpg_utils.config import get_config
from cpg_pipes.pipeline import pipeline_entry_point
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.genotype_sample import GenotypeSample
from cpg_pipes.utils import exists

logger = logging.getLogger(__file__)


@pipeline_entry_point(name='sample genotyping')
def main(pipeline: Pipeline):
    """
    Batch pipeline to generate GVCF for all samples CRAMs.
    """
    if get_config()['workflow'].get('skip_samples_with_missing_input'):
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(f'Could not find CRAM, skipping sample {sample.id}')
                sample.active = False

    pipeline.run(
        stages=[GenotypeSample],
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
