"""
Stage to run STR analysis with STRipy-pipeline.

See https://gitlab.com/andreassh/stripy-pipeline
"""

from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(analysis_type='other')
class DevStage1(SequencingGroupStage):
    """
    A super fun dev stage
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        j = get_batch().new_job('devstage1')
        j.command(f'echo "{sequencing_group.id} says hello from devstage1"')
        jobs.append(j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
