"""
Test building stages DAG.
"""

import pytest
from pytest_mock import MockFixture

from .. import set_config
from . import run_workflow


def test_first_last_stages_misconfigured(mocker: MockFixture, tmp_path):
    """
    A -> B -> C -> D
    first_stages = [C]
    last_stages = [B]
    Will raise no stages to run
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

    first_stages = ['C']
    last_stages = ['B']

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
    from cpg_workflows.workflow import WorkflowError

    with pytest.raises(WorkflowError, match='No stages to run'):
        run_workflow(mocker)
