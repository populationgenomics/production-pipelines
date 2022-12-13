"""
Test workflow status reporter.
"""

import toml
import pytest
from pytest_mock import MockFixture
from cpg_utils import to_path, Path
from . import results_prefix, update_dict

TOML = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = '<stub>'
sequencing_type = 'genome'
status_reporter = 'metamist'

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

[storage.default]
default = '{results_prefix()}'

[storage.fewgenomes]
default = '{results_prefix()}'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
dry_run = true
"""


def _mock_config() -> dict:
    d: dict = {}
    for fp in [
        to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml',
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'seqr_loader.toml',
    ]:
        with fp.open():
            update_dict(d, toml.load(fp))

    update_dict(d, toml.loads(TOML))
    return d


def _common(mocker):
    mocker.patch('cpg_utils.config.get_config', _mock_config)

    def mock_create_new_analysis(_, project, analysis_model) -> int:
        print(f'Analysis model in project {project}: {analysis_model}')
        return 1  # metamist "analysis" entry ID

    mocker.patch(
        'sample_metadata.apis.AnalysisApi.create_new_analysis', mock_create_new_analysis
    )

    from cpg_workflows.targets import Cohort

    def mock_create_cohort() -> Cohort:
        c = Cohort()
        ds = c.create_dataset('my_dataset')
        ds.add_sample('CPG01', external_id='SAMPLE1')
        ds.add_sample('CPG02', external_id='SAMPLE2')
        return c

    mocker.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)


def test_status_reporter(mocker: MockFixture):
    _common(mocker)

    from cpg_utils.hail_batch import dataset_path
    from cpg_workflows.inputs import get_cohort
    from cpg_workflows.targets import Sample
    from cpg_workflows.batch import get_batch
    from cpg_workflows.workflow import (
        SampleStage,
        StageInput,
        StageOutput,
        stage,
        run_workflow,
    )

    @stage(analysis_type='qc')
    class MyQcStage1(SampleStage):
        """
        Just a sample-level stage.
        """

        def expected_outputs(self, sample: Sample) -> Path:
            return to_path(dataset_path(f'{sample.id}.tsv'))

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = get_batch().new_job(
                'Echo', self.get_job_attrs(sample) | dict(tool='echo')
            )
            j.command(f'echo {sample.id}_done >> {j.output}')
            get_batch().write_output(j.output, str(self.expected_outputs(sample)))
            print(f'Writing to {self.expected_outputs(sample)}')
            return self.make_outputs(sample, self.expected_outputs(sample), [j])

    @stage(analysis_type='qc', analysis_key='bed')
    class MyQcStage2(SampleStage):
        """
        Just a sample-level stage.
        """

        def expected_outputs(self, sample: Sample) -> dict:
            return {
                'bed': to_path(dataset_path(f'{sample.id}.bed')),
                'tsv': to_path(dataset_path(f'{sample.id}.tsv')),
            }

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = get_batch().new_job(
                'Echo', self.get_job_attrs(sample) | dict(tool='echo')
            )
            j.command(f'echo {sample.id}_done >> {j.output}')
            get_batch().write_output(
                j.output, str(self.expected_outputs(sample)['bed'])
            )
            print(f'Writing to {self.expected_outputs(sample)["bed"]}')
            return self.make_outputs(sample, self.expected_outputs(sample), [j])

    run_workflow(stages=[MyQcStage1, MyQcStage2])

    assert 'metamist' in get_batch().job_by_tool, get_batch().job_by_tool
    assert (
        get_batch().job_by_tool['metamist']['job_n']
        == len(get_cohort().get_samples()) * 4
    )


def test_status_reporter_fails(mocker: MockFixture):
    _common(mocker)

    from cpg_utils.hail_batch import dataset_path
    from cpg_workflows.targets import Sample
    from cpg_workflows.batch import get_batch
    from cpg_workflows.workflow import (
        SampleStage,
        StageInput,
        StageOutput,
        stage,
        run_workflow,
    )

    @stage(analysis_type='qc')
    class MyQcStage(SampleStage):
        """
        Just a sample-level stage.
        """

        def expected_outputs(self, sample: Sample) -> dict:
            return {
                'bed': dataset_path(f'{sample.id}.bed'),
                'tsv': dataset_path(f'{sample.id}.tsv'),
            }

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = get_batch().new_job(
                'Echo', self.get_job_attrs(sample) | dict(tool='echo')
            )
            j.command(f'echo {sample.id}_done >> {j.output}')
            get_batch().write_output(
                j.output, str(self.expected_outputs(sample)['bed'])
            )
            print(f'Writing to {self.expected_outputs(sample)["bed"]}')
            return self.make_outputs(sample, self.expected_outputs(sample), [j])

    from cpg_workflows.workflow import WorkflowError

    with pytest.raises(WorkflowError):
        run_workflow(stages=[MyQcStage])
