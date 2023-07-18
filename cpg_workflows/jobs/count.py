"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)


class FeatureCounts:
    """
    Construct a featureCounts command for counting reads.
    """

    def __init__(self):
        self.command = ['featureCounts']

    def __str__(self):
        return ' '.join(self.command)
    
    def __repr__(self):
        return self.__str__()


def count(
    b: hb.Batch,
    input_bam: BamPath,
    output_path: str | Path,
):
    pass
