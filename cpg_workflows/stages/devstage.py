"""
Example Stages
"""

from dataclasses import dataclass

from cpg_utils.config import config_retrieve
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

    # put the images here, so you're not calling image_path later on
    # everything you need for the stage (except stage input) should be in this config
    image_samtools: str

    reference_version: str | None = None
    resource_override: str | None = None
    skip_stage: bool = False


class ExampleStageOutputs(StageOutputs):
    test_output: str


@stage(analysis_type='custom', analysis_keys=['test_output'])
class ExampleStage(SequencingGroupStage):
    """
    An example stage
    """

    def setup_config_inputs(self, global_config: dict) -> ExampleStageConfig:
        return ExampleStageConfig(
            reference_version=config_retrieve(['workflow', 'reference_version'], global_config),
            images_samtools=config_retrieve(['workflow', 'images', 'samtools'], global_config),
        )

    def expected_outputs(self, sequencing_group: SequencingGroup) -> ExampleStageOutputs:
        return ExampleStageOutputs(
            test_output=sequencing_group.dataset.prefix() / 'example' / f'{sequencing_group.id}_example.txt',
        )

    def queue_jobs(
        self,
        sequencing_group: SequencingGroup,
        inputs: StageInput,
        stage_config: ExampleStageConfig,
    ) -> StageOutput | None:
        jobs = []

        # no get_batch(), batch ref is provided as instance on class
        j = example_job.example_job(
            batch=self.batch,
            sequencing_group=sequencing_group,
            config=stage_config,
        )
        jobs.append(j)
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
