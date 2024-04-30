"""
Test building stages DAG.
"""

from pytest_mock import MockFixture

from .. import set_config
from . import run_workflow


def test_force_stages(mocker: MockFixture, tmp_path):
    """
    A -> B -> C, all results exist
    force_stages = [B2]
    Should run: B
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

    # Skip stages with outputs that already exist
    check_expected_outputs = true

    # Except force stage B to re-run
    force_stages = ['B']

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

    # dummy mocking to avoid file system scanning
    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: True)
    set_config(conf, tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch

    run_workflow(mocker)
    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B']['job_n'] == 1
    assert 'C' not in get_batch().job_by_stage
