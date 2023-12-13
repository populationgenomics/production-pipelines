"""
Test building stages DAG.
"""

from pytest_mock import MockFixture

from .. import set_config
from . import run_workflow


def test_first_last_stages(mocker: MockFixture, tmp_path):
    """
    A -> B -> C -> D
    first_stages = [B]
    last_stages = [C]
    Should run: B -> C
    """
    conf = f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    driver_image = 'stub'

    check_inputs = false
    check_intermediates = false
    check_expected_outputs = false

    first_stages = ['B']
    last_stages = ['C']

    [storage.default]
    default = '{tmp_path}'
    [storage.fewgenomes]
    default = '{tmp_path}'

    [hail]
    billing_project = 'fewgenomes'
    delete_scratch_on_exit = false
    backend = 'local'
    dry_run = true
    """
    set_config(conf, tmp_path / 'config.toml')
    run_workflow(mocker)

    from cpg_utils.hail_batch import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B']['job_n'] == 1
    assert get_batch().job_by_stage['C']['job_n'] == 1
    assert 'D' not in get_batch().job_by_stage
