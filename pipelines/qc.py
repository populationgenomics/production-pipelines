#!/usr/bin/env python3

"""
Batch pipeline to run QC on all CRAMs. Outputs MultiQC reports for each dataset, 
and pedigree checks.
"""

import logging

from cpg_pipes.pipeline import pipeline_entry_point
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.multiqc import MultiQC
from cpg_pipes.utils import exists
from cpg_utils.config import get_config

logger = logging.getLogger(__file__)


@pipeline_entry_point(name='CRAM QC')
def main(pipeline: Pipeline):
    """
    Batch pipeline to run QC on all CRAMs. Outputs MultiQC reports for each dataset, 
    and pedigree checks.
    """
    if get_config()['workflow'].get('skip_samples_with_missing_input'):
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(
                    f'Could not find CRAM, skipping sample {sample.id}'
                )
                sample.active = False

    pipeline.run(
        stages=[MultiQC],
        force_all_implicit_stages=True,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
