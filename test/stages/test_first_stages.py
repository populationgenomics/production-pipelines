"""
Test building stages DAG.
"""

from pytest_mock import MockFixture

from .. import set_config
from . import run_workflow


def test_first_stages(mocker: MockFixture, tmp_path):
    """
    A -> B -> C
    A2 -> B2 -> C2
    first_stages = [B2]
    Should run: B2 -> C2
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
    check_expected_outputs = true

    first_stages = ['B2']

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
    from cpg_utils.hail_batch import get_batch

    run_workflow(mocker)
    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert 'B' not in get_batch().job_by_stage
    assert 'C' not in get_batch().job_by_stage
    assert 'A2' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B2']['job_n'] == 1
    assert get_batch().job_by_stage['C2']['job_n'] == 1
