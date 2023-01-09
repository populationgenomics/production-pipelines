"""
Test initialising Batch object.
"""

import toml
from pytest_mock import MockFixture
from cpg_utils import to_path
from . import results_prefix

TOML = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
sequencing_type = 'genome'

check_inputs = false
check_intermediates = false
check_expected_outputs = false

[storage.default]
default = '{results_prefix()}'

[storage.fewgenomes]
default = '{results_prefix()}'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
"""


def test_batch_job(mocker: MockFixture):
    """
    Test creating a job and running a batch.
    """
    mocker.patch('cpg_utils.config.get_config', lambda: toml.loads(TOML))

    from cpg_utils.hail_batch import dataset_path, command
    from cpg_workflows.batch import get_batch

    b = get_batch('Test batch job')
    j1 = b.new_job('Jo b1')
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
