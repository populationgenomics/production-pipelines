"""
Example Stages
"""

from dataclasses import dataclass

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import example_job
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(analysis_type='custom', analysis_keys=['test'])
class DevStage1(SequencingGroupStage):
    """
    A super fun dev stage
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {'test': sequencing_group.dataset.prefix() / 'dev' / f'{sequencing_group.id}_devstage1.txt'}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        j = get_batch().new_job('devstage1')
        j.command(f'echo "{sequencing_group.id} says hello from devstage1" > {j.out}')
        jobs.append(j)
        get_batch().write_output(
            j.out,
            str(sequencing_group.dataset.prefix() / 'dev' / f'{sequencing_group.id}_devstage1.txt'),
        )

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@stage(analysis_type='custom', analysis_keys=['test_output'])
class ExampleStagex(SequencingGroupStage):
    """
    An example stage
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {'test_output': sequencing_group.dataset.prefix() / 'example' / f'{sequencing_group.id}_example.txt'}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        j = example_job.example_job(
            batch=get_batch(),
            sequencing_group=sequencing_group,
        )
        jobs.append(j)
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@dataclass
class ExampleStageConfig(StageConfig):
    """
    An example stage config
    """

    reference_version: str | None = None
    resource_override: str | None = None
    skip_stage: bool = False


@stage(analysis_type='custom', analysis_keys=['test_output'])
class ExampleStage(SequencingGroupStage):
    """
    An example stage
    """

    def setup_config_inputs(example_stage_config: dict) -> ExampleStageConfig:
        return ExampleStageConfig(**example_stage_config)

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {'test_output': sequencing_group.dataset.prefix() / 'example' / f'{sequencing_group.id}_example.txt'}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        stage_config = self.setup_config_inputs(get_config('example_stage'))
        j = example_job.example_job(batch=get_batch(), sequencing_group=sequencing_group, config=stage_config)
        jobs.append(j)
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
