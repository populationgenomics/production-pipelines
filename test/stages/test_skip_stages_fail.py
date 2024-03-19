"""
Test building stages DAG.
"""

import pytest
from pytest_mock import MockFixture

from .. import set_config
from . import run_workflow


def test_skip_stages_fail(mocker: MockFixture, tmp_path):
    """
    A -> B -> C
    A2 -> B2 -> C2
    skip_stages = [A2]
    check_expected_outputs = True
    Should raise WorkflowError (A2 outputs don't exist)
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

    skip_stages = ['A2']

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

    with pytest.raises(
        WorkflowError,
        match='A2: stage is required, but is skipped, and the following expected outputs',
    ):
        run_workflow(mocker)
