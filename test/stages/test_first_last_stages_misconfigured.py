"""
Test building stages DAG.
"""
import pytest
import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_first_last_stages_misconfigured(mocker: MockFixture):
    """
    A -> B -> C -> D
    first_stages = [C]
    last_stages = [B]
    Will raise no stages to run
    """
    conf = toml.loads(TOML)
    conf['workflow']['first_stages'] = ['C']
    conf['workflow']['last_stages'] = ['B']
    conf['workflow']['check_expected_outputs'] = False
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    from cpg_workflows.workflow import WorkflowError

    with pytest.raises(
        WorkflowError,
        match='No stages to run',
    ):
        run_workflow(mocker)
