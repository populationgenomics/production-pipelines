"""
Test initialising Batch object.
"""

from cpg_utils import to_path

from . import set_config

from pathlib import Path
from cpg_utils import to_path
from pytest_mock import MockFixture

from cpg_utils.hail_batch import dataset_path

from cpg_workflows.batch import get_batch
from cpg_workflows.inputs import get_cohort
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    run_workflow,
    stage,
)
from cpg_workflows.targets import Cohort

CONFIG = """
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
default = '{directory}'

[storage.fewgenomes]
default = '{directory}'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
dry_run = true

[images]
cpg_workflows = "stub"
"""


def test_batch_job(tmp_path):
    """
    Test creating a job and running a batch.
    """
    config = f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    driver_image = 'test'

    check_inputs = false
    check_intermediates = false
    check_expected_outputs = false

    [storage.default]
    default = '{tmp_path}'

    [storage.fewgenomes]
    default = '{tmp_path}'

    [hail]
    billing_project = 'fewgenomes'
    delete_scratch_on_exit = false
    backend = 'local'
    """
    set_config(config, tmp_path / 'config.toml')

    from cpg_utils.hail_batch import command, dataset_path

    from cpg_workflows.batch import get_batch

    b = get_batch('Test batch job')
    j1 = b.new_job('Job 1')
    text = 'success'
    cmd = f"""\
    echo {text} > {j1.output}
    """
    j1.command(command(cmd))
    output1_path = dataset_path('output1.txt')
    b.write_output(j1.output, str(output1_path))

    j2 = b.new_job('Job 2')
    j2.command(f'touch {j2.output}')
    j2.command(f'cat {b.read_input(output1_path)} >> {j2.output}')
    j2.depends_on(j1)
    output2_path = dataset_path('output2.txt')
    b.write_output(j2.output, str(output2_path))

    b.run()
    with to_path(output2_path).open() as fh:
        assert fh.read().strip() == text


def mock_create_analysis(_, project, analysis_model) -> int:
    print(f'Analysis model in project {project}: {analysis_model}')
    return 1  # metamist "analysis" entry ID


def mock_create_cohort() -> Cohort:
    c = Cohort()
    ds = c.create_dataset('my_dataset')
    ds.add_sequencing_group('CPG01', external_id='SAMPLE1')
    ds.add_sequencing_group('CPG02', external_id='SAMPLE2')
    return c


@stage(analysis_type='qc')
class MyQcStage1(SequencingGroupStage):
    """
    A dummy SequencingGroup Stage.
    """

    def expected_outputs(self, sample: SequencingGroup) -> Path:
        return to_path(dataset_path(f'{sample.id}.tsv'))

    def queue_jobs(
        self, sample: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        j = get_batch().new_job('Echo', self.get_job_attrs(sample) | dict(tool='echo'))
        j.command(f'echo {sample.id}_done >> {j.output}')
        get_batch().write_output(j.output, str(self.expected_outputs(sample)))
        print(f'Writing to {self.expected_outputs(sample)}')
        return self.make_outputs(sample, self.expected_outputs(sample), [j])


@stage(analysis_type='qc', analysis_keys=['bed'])
class MyQcStage2(SequencingGroupStage):
    """
    A dummy SequencingGroup stage.
    """

    def expected_outputs(self, sample: SequencingGroup) -> dict:
        return {
            'bed': to_path(dataset_path(f'{sample.id}.bed')),
            'tsv': to_path(dataset_path(f'{sample.id}.tsv')),
        }

    def queue_jobs(
        self, sample: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        j = get_batch().new_job('Echo', self.get_job_attrs(sample) | dict(tool='echo'))
        j.command(f'echo {sample.id}_done >> {j.output}')
        get_batch().write_output(j.output, str(self.expected_outputs(sample)['bed']))
        print(f'Writing to {self.expected_outputs(sample)["bed"]}')
        return self.make_outputs(sample, self.expected_outputs(sample), [j])


def test_attributes(mocker: MockFixture, tmp_path):
    mocker.patch('metamist.apis.AnalysisApi.create_analysis', mock_create_analysis)
    mocker.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)

    set_config(CONFIG, tmp_path / 'config.toml')

    workflow_stages = [MyQcStage1, MyQcStage2]
    run_workflow(stages=workflow_stages)

    # Check that the correct number of jobs were created
    assert (
        len(get_batch()._jobs)
        == len(get_cohort().get_sequencing_groups()) * len(workflow_stages) * 3
    )
    # 3 jobs per stage, assumes no nested workflow stages

    # Check that all expected attributes are present
    expected_attrs = [
        'stage',
        'sequencing_type',
        'dataset',
        'sequencing_group',
        'participant_id',
        'tool',
        'sequencing_groups',
    ]

    for job in get_batch()._jobs:
        for key in job.attributes:
            assert key in expected_attrs
        assert job.attributes['stage'] in [s.__name__ for s in workflow_stages]
        assert job.attributes['sequencing_type'] == 'genome'
        assert job.attributes['dataset'] == 'my_dataset'
        assert job.attributes['tool'] in ['echo', 'metamist']
        assert job.attributes['participant_id'] in ['SAMPLE1', 'SAMPLE2']
        assert job.attributes['sequencing_group'] in ['CPG01', 'CPG02']
        assert (
            job.attributes['sequencing_groups'] == "['CPG01']"
            or job.attributes['sequencing_groups'] == "['CPG02']"
        )

    # Check that the job_by_stage and job_by_tool dicts are correct
    for stage, job in get_batch().job_by_stage.items():
        assert stage in [s.__name__ for s in workflow_stages]
        assert job['sequencing_groups'] == {'CPG01', 'CPG02'}

    for tool, job in get_batch().job_by_tool.items():
        assert tool in ['echo', 'metamist']
        assert job['sequencing_groups'] == {'CPG01', 'CPG02'}
