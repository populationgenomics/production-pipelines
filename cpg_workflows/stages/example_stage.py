"""
Example Stages
"""

from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, get_workflow, stage


@stage(analysis_type='custom', analysis_keys=['test'])
class ExampleStage1(SequencingGroupStage):
    """
    A super fun dev stage
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {
            'test': sequencing_group.dataset.prefix()
            / 'dev'
            / get_workflow().output_version
            / f'{sequencing_group.id}_devstage1.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        j = get_batch().new_job('examplestage1')
        j.command(f'echo "{sequencing_group.id} says hello from example_stage_1" > {j.out}')
        jobs.append(j)
        get_batch().write_output(
            j.out,
            str(
                sequencing_group.dataset.prefix()
                / 'dev'
                / get_workflow().output_version
                / f'{sequencing_group.id}_devstage1.txt',
            ),
        )

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@stage(required_stages=[ExampleStage1], analysis_type='custom', analysis_keys=['test_output'])
class ExampleStage2(SequencingGroupStage):
    """
    An example stage
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {
            'test_output': sequencing_group.dataset.prefix()
            / 'example'
            / get_workflow().output_version
            / f'{sequencing_group.id}_example.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = []
        j = get_batch().new_job('examplestage1')
        j.command(f'echo "{sequencing_group.id} says hello from example_stage_2" > {j.out}')
        jobs.append(j)
        get_batch().write_output(
            j.out,
            str(
                sequencing_group.dataset.prefix()
                / 'example'
                / get_workflow().output_version
                / f'{sequencing_group.id}_example.txt',
            ),
        )

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
