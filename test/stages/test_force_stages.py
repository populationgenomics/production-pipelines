"""
Test building stages DAG.
"""

import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_force_stages(mocker: MockFixture):
    """
    A -> B -> C, all results exist
    force_stages = [B2]
    Should run: B
    """
    conf = toml.loads(TOML)
    conf['workflow']['force_stages'] = ['B']
    mocker.patch('cpg_utils.config.get_config', lambda: conf)
    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: True)

    run_workflow(mocker)

    from cpg_workflows.workflow import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B']['job_n'] == 1
    assert 'C' not in get_batch().job_by_stage
