"""
Test building stages DAG.
"""

import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_first_stages(mocker: MockFixture):
    """
    A -> B -> C
    A2 -> B2 -> C2
    first_stages = [B2]
    Should run: B2 -> C2
    """
    conf = toml.loads(TOML)
    conf['workflow']['first_stages'] = ['B2']
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    run_workflow(mocker)

    from cpg_workflows.workflow import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert 'B' not in get_batch().job_by_stage
    assert 'C' not in get_batch().job_by_stage
    assert 'A2' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B2']['job_n'] == 1
    assert get_batch().job_by_stage['C2']['job_n'] == 1
