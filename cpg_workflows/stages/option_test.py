"""
This is a test to see if we can create a pipeline with optional stages
"""

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage

STAGE_OPTIONS = ['continue', 'replace']

@stage()
class DataSource(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'data': self.tmp_prefix / 'data.txt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        data_source_job = get_batch().new_job('DataSource', self.get_job_attrs())
        data_source_job.image(config_retrieve(['workflow', 'driver_image']))
        data_source_job.command(f'echo "Hello, World!" > {data_source_job.output}')
        get_batch().write_output(data_source_job.output, str(outputs['data']))

        return self.make_outputs(multicohort, data=outputs, jobs=data_source_job)


@stage(required_stages=[DataSource])
class DataContinue(MultiCohortStage):
    """
    this stage takes the original data and passes it on
    """
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'data': self.tmp_prefix / 'continue.txt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(multicohort)

        if config_retrieve(['workflow', 'second_stage']) != 'continue':
            return self.make_outputs(multicohort, data=outputs)

        input_data = get_batch().read_input(str(inputs.as_path(target=multicohort, stage=DataSource, key='data')))

        data_continue_job = get_batch().new_job('DataContinue', self.get_job_attrs())
        data_continue_job.image(config_retrieve(['workflow', 'driver_image']))

        # take the input, and write it to the output
        data_continue_job.command(f'mv {input_data} {data_continue_job.output}')
        get_batch().write_output(data_continue_job.output, str(outputs['data']))
        return self.make_outputs(multicohort, data=outputs, jobs=data_continue_job)


@stage(required_stages=[DataSource])
class DataReplace(MultiCohortStage):
    """
    this stage replaces that data with something different
    """
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'data': self.tmp_prefix / 'replace.txt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(multicohort)

        if config_retrieve(['workflow', 'second_stage']) != 'replace':
            return self.make_outputs(multicohort, data=outputs)

        data_replace_job = get_batch().new_job('DataReplace', self.get_job_attrs())
        data_replace_job.image(config_retrieve(['workflow', 'driver_image']))
        data_replace_job.command(f'echo "Goodbye, Cruel World!" > {data_replace_job.output}')
        get_batch().write_output(data_replace_job.output, str(outputs['data']))
        return self.make_outputs(multicohort, data=outputs, jobs=data_replace_job)


@stage(required_stages=[DataReplace, DataContinue])
class DataSink(MultiCohortStage):
    """
    this stage takes the data from one of the two previous stages and prints
    """
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        file_name = config_retrieve(['workflow', 'second_stage']) + '.txt'
        return {'data': self.tmp_prefix / file_name}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        stage_to_use = config_retrieve(['workflow', 'second_stage'])

        if stage_to_use == 'continue':
            input_data = get_batch().read_input(str(inputs.as_path(target=multicohort, stage=DataContinue, key='data')))
        elif stage_to_use == 'replace':
            input_data = get_batch().read_input(str(inputs.as_path(target=multicohort, stage=DataReplace, key='data')))
        else:
            raise ValueError(f'config.workflow.second_stage must be one of {STAGE_OPTIONS}')

        outputs = self.expected_outputs(multicohort)

        final_job = get_batch().new_job('DataSink', self.get_job_attrs())
        final_job.command(f'cat {input_data}')
        final_job.command(f'mv {input_data} {final_job.output}')

        get_batch().write_output(final_job.output, str(outputs['data']))

        return self.make_outputs(multicohort, data=outputs, jobs=final_job)
