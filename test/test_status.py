"""
Test workflow status reporter.
"""

import toml
from pytest_mock import MockFixture

from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths, update_dict
from cpg_utils.hail_batch import dataset_path
from cpg_workflows.inputs import get_cohort
from cpg_workflows.targets import Sample, Cohort
from cpg_workflows.utils import timestamp
from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import (
    SampleStage,
    StageInput,
    StageOutput,
    ExpectedResultT,
    stage,
    run_workflow,
)

tmp_dir_path = to_path(__file__).parent / 'results' / timestamp()
tmp_dir_path = tmp_dir_path.absolute()
tmp_dir_path.mkdir(parents=True, exist_ok=True)

DEFAULT_CONF = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = '<stub>'
sequencing_type = 'genome'

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
"""


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    d['workflow']['local_dir'] = str(dir_path)
    if extra_conf:
        update_dict(d, extra_conf)
    conf_path = dir_path / 'config.toml'
    with conf_path.open('w') as f:
        toml.dump(d, f)

    set_config_paths(
        [
            str(p)
            for p in [
                to_path(__file__).parent.parent
                / 'configs'
                / 'defaults'
                / 'workflows.toml',
                to_path(__file__).parent.parent
                / 'configs'
                / 'defaults'
                / 'seqr_loader.toml',
                conf_path,
            ]
        ]
    )


def test_status_reporter(mocker: MockFixture):
    """
    Testing metamist status reporter.
    """
    _set_config(
        tmp_dir_path,
        {
            'workflow': {
                'status_reporter': 'metamist',
            },
            'hail': {
                'dry_run': True,
            },
        },
    )

    def mock_create_new_analysis(_, project, analysis_model) -> int:
        print(f'Analysis model in project {project}: {analysis_model}')
        return 1  # metamist "analysis" entry ID

    mocker.patch(
        'sample_metadata.apis.AnalysisApi.create_new_analysis', mock_create_new_analysis
    )

    def mock_create_cohort() -> Cohort:
        c = Cohort()
        ds = c.create_dataset('my_dataset')
        ds.add_sample('CPG01', external_id='SAMPLE1')
        ds.add_sample('CPG02', external_id='SAMPLE2')
        return c

    mocker.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)

    @stage(analysis_type='qc')
    class MyQcStage(SampleStage):
        """
        Just a sample-level stage.
        """

        def expected_outputs(self, sample: Sample) -> ExpectedResultT:
            return dataset_path(f'{sample.id}.tsv')

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = self.b.new_job('Echo', self.get_job_attrs(sample) | dict(tool='echo'))
            j.command(f'echo {sample.id}_done >> {j.output}')
            self.b.write_output(j.output, str(self.expected_outputs(sample)))
            print(f'Writing to {self.expected_outputs(sample)}')
            return self.make_outputs(sample, self.expected_outputs(sample), [j])

    run_workflow(stages=[MyQcStage])

    assert 'metamist' in get_batch().job_by_tool, get_batch().job_by_tool
    assert (
        get_batch().job_by_tool['metamist']['job_n']
        == len(get_cohort().get_samples()) * 2
    )
