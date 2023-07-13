"""
concept test of a stage using a test_mode
"""
import inspect
from cpg_workflows.stages.happy_validation import ValidationMtToVcf
from cpg_workflows.query_modules.validation import single_sample_vcf_from_dataset_vcf
from unittest import mock

from cpg_workflows.batch import get_batch

from .. import set_config

from test.test_workflow import mock_create_cohort


TOML = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = 'test'
sequencing_type = 'genome'

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

[storage.default]
default = "{directory}"

[storage.fewgenomes]
default = "{directory}"

[storage.my_dataset]
default = "{directory}"

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'

[references.SAMPLE1]
[references.SAMPLE2]
[inputs]
sample_hash = "boop"

[images]
cpg_workflows = "test"
"""


@mock.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)
@mock.patch('cpg_workflows.workflow.list_of_all_dir_contents', set)
def test_validation_mt_to_vcf(tmp_path):

    from cpg_workflows.workflow import run_workflow

    conf = TOML.format(directory=tmp_path)
    set_config(conf, tmp_path / 'config.toml')
    run_workflow(stages=[ValidationMtToVcf], test_mode=True)
    jobs = get_batch()._jobs

    expected_code = inspect.getsource(single_sample_vcf_from_dataset_vcf)
    assert len(jobs) == 2
    for job in jobs:
        assert job.attributes['stage'] == 'ValidationMtToVcf'
        assert expected_code in job._command[0]
