"""
Test building Workflow object.
"""

import toml
from pytest_mock import MockFixture

from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths, update_dict
from cpg_utils.hail_batch import dataset_path
from cpg_workflows import get_batch
from cpg_workflows.inputs import get_cohort
from cpg_workflows.targets import Sample, Cohort
from cpg_workflows.utils import timestamp
from cpg_workflows.workflow import (
    SampleStage,
    StageInput,
    StageOutput,
    CohortStage,
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
    config_path = dir_path / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(config_path)])


def test_workflow(mocker: MockFixture):
    """
    Testing running a workflow from a mock cohort.
    """
    _set_config(tmp_dir_path)

    def mock_create_cohort() -> Cohort:
        c = Cohort()
        ds = c.create_dataset('my_dataset')
        ds.add_sample('CPG01', external_id='SAMPLE1')
        ds.add_sample('CPG02', external_id='SAMPLE2')
        return c

    mocker.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)

    output_path = to_path(dataset_path('cohort.tsv'))

    assert len(get_cohort().get_samples()) == 2

    @stage
    class MySampleStage(SampleStage):
        """
        Just a sample-level stage.
        """

        def expected_outputs(self, sample: Sample) -> Path:
            return to_path(dataset_path(f'{sample.id}.tsv'))

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = get_batch().new_job('Sample job', self.get_job_attrs(sample))
            j.command(f'echo {sample.id}_done >> {j.output}')
            get_batch().write_output(j.output, str(self.expected_outputs(sample)))
            print(f'Writing to {self.expected_outputs(sample)}')
            return self.make_outputs(sample, self.expected_outputs(sample))

    @stage(required_stages=MySampleStage)
    class MyCohortStage(CohortStage):
        """
        Just a cohort-level stage.
        """

        def expected_outputs(self, cohort: Cohort) -> Path:
            return output_path

        def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
            path_by_sample = inputs.as_path_by_target(MySampleStage)
            assert len(path_by_sample) == len(cohort.get_samples())
            j = get_batch().new_job('Cohort job', self.get_job_attrs(cohort))
            j.command(f'touch {j.output}')
            for _, sample_result_path in path_by_sample.items():
                input_file = get_batch().read_input(str(sample_result_path))
                j.command(f'cat {input_file} >> {j.output}')
            get_batch().write_output(j.output, str(self.expected_outputs(cohort)))
            print(f'Writing to {self.expected_outputs(cohort)}')
            return self.make_outputs(cohort, self.expected_outputs(cohort))

    run_workflow(stages=[MyCohortStage])

    print(f'Checking result in {output_path}:')
    with output_path.open() as f:
        result = f.read()
        print(result)
        assert result.split() == ['CPG01_done', 'CPG02_done'], result
