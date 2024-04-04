"""
Stage to run STR analysis with STRipy-pipeline.

See https://gitlab.com/andreassh/stripy-pipeline
"""

from cpg_utils import Path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)

from cpg_workflows.jobs import first_class


@stage(analysis_keys=['first_class_out'], analysis_type='first_class')
# @stage()
class FirstClass(SequencingGroupStage):
    """
    Test First Class Files
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:

        return {
            'first_class_out': f'output/{sequencing_group.id}/first_class.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:

        jobs = []

        j = first_class.first_class(
            b=get_batch(),
            sequencing_group=sequencing_group,
            out_path=self.expected_outputs(sequencing_group)['first_class_out'],
        )
        jobs.append(j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
