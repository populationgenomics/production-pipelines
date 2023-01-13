"""
Test building stages DAG.
"""
import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_skip_stages(mocker: MockFixture):
    """
    A -> B -> C
    A2 -> B2 -> C2
    skip_stages = [A2]
    Should run: A -> B -> C, B2 -> C2
    """
    conf = toml.loads(TOML)
    conf['workflow']['skip_stages'] = ['A2']
    conf['workflow']['check_expected_outputs'] = False
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    run_workflow(mocker)

    from cpg_workflows.workflow import get_batch

    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage['A']['job_n'] == 1
    assert get_batch().job_by_stage['B']['job_n'] == 1
    assert get_batch().job_by_stage['C']['job_n'] == 1
    assert 'A2' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B2']['job_n'] == 1
    assert get_batch().job_by_stage['C2']['job_n'] == 1
