"""
Perform outlier gene expression analysis with Outrider.
"""

import logging
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
)
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.stages.count import Count
from cpg_workflows.jobs import outrider


@stage(
    required_stages=Count,
)
class Outrider(SequencingGroupStage):
    """
    Perform outlier gene expression analysis with Outrider.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate outrider outputs.
        """
        return {
            'test': sequencing_group.dataset.prefix() / 'outrider' / f'{sequencing_group.id}.test',
        }
    
    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run outrider.
        """
        j = outrider.outrider(
            b=get_batch(),
            output_path=self.expected_outputs(sequencing_group)['test'],
            sample_name=sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=j)
