"""
Test building stages DAG.
"""

import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_last_stages(mocker: MockFixture):
    """
    A -> B -> C
    A2 -> B2 -> C2
    last_stages = [B2]
    Should run: A2 -> B2
    """
    conf = toml.loads(TOML)
    conf['workflow']['last_stages'] = ['B2']
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    run_workflow(mocker)

    from cpg_workflows.workflow import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert 'B' not in get_batch().job_by_stage
    assert 'C' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['A2']['job_n'] == 1
    assert get_batch().job_by_stage['B2']['job_n'] == 1
    assert 'C2' not in get_batch().job_by_stage
