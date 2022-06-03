#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

import logging

from cpg_utils.config import get_config

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.utils import exists
from cpg_pipes.pipeline import pipeline_entry_point
from cpg_pipes.stages.somalier import CramSomalierAncestry, CramSomalierPedigree

logger = logging.getLogger(__file__)


@pipeline_entry_point(name='pedigree')
def main(pipeline: Pipeline):
    """
    Entry point, decorated by pipeline click options.
    """
    if get_config()['workflow'].get('skip_samples_with_missing_input'):
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(f'Could not find CRAM, skipping sample {sample.id}')
                sample.active = False

    pipeline.run(
        stages=[
            CramSomalierPedigree,
            CramSomalierAncestry,
        ],
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
