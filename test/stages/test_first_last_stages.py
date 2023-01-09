"""
Test building stages DAG.
"""

import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_first_last_stages(mocker: MockFixture):
    """
    A -> B -> C -> D
    first_stages = [B]
    last_stages = [C]
    Should run: B -> C
    """
    conf = toml.loads(TOML)
    conf['workflow']['first_stages'] = ['B']
    conf['workflow']['last_stages'] = ['C']
    conf['workflow']['check_expected_outputs'] = False
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    run_workflow(mocker)

    from cpg_workflows.workflow import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B']['job_n'] == 1
    assert get_batch().job_by_stage['C']['job_n'] == 1
    assert 'D' not in get_batch().job_by_stage
