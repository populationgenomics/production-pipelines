"""
Test building stages DAG.
"""
import pytest
import toml
from pytest_mock import MockFixture
from . import TOML, run_workflow


def test_skip_stages_fail(mocker: MockFixture):
    """
    A -> B -> C
    A2 -> B2 -> C2
    skip_stages = [A2]
    check_expected_outputs = True
    Should raise WorkflowError (A2 outputs don't exist)
    """
    conf = toml.loads(TOML)
    conf['workflow']['skip_stages'] = ['A2']
    conf['workflow']['check_expected_outputs'] = True
    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    from cpg_workflows.workflow import WorkflowError

    with pytest.raises(
        WorkflowError,
        match='A2: stage is required, but is skipped, and the following expected outputs',
    ):
        run_workflow(mocker)
