"""
Test building stages DAG.
"""

import pytest
from pytest_mock import MockFixture

from cpg_utils.hail_batch import get_batch, reset_batch
from cpg_workflows.workflow import WorkflowError, stage

from .. import set_config
from . import TestStage, run_workflow


@stage
class A(TestStage):
    pass


@stage(required_stages=A)
class B1(TestStage):
    pass


@stage(required_stages=A)
class B2(TestStage):
    pass


@stage(required_stages=B1)
class C1(TestStage):
    pass


@stage(required_stages=B2)
class C2(TestStage):
    pass


def create_config(tmp_path, allow_missing_outputs_for_stages=None):
    return f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    driver_image = 'stub'

    check_inputs = false
    check_intermediates = false
    check_expected_outputs = true

    allow_missing_outputs_for_stages = {allow_missing_outputs_for_stages or []}
    only_stages = ['B1']

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


def test_works_when_using_allow_missing_outputs_for_stages_on_required_parent_stage(mocker: MockFixture, tmp_path):
    """
    A -> B1 -> C1
    A -> B2 -> C2

    only_stages = [B1]
    Should run: B1
    """
    conf = create_config(tmp_path, allow_missing_outputs_for_stages=['A'])
    set_config(conf, tmp_path / 'config.toml')
    run_workflow(mocker, stages=[C1, C2])
    print('Job by stage:', get_batch().job_by_stage)
    assert 'A' not in get_batch().job_by_stage
    assert get_batch().job_by_stage['B1']['job_n'] == 1
    assert 'B2' not in get_batch().job_by_stage
    assert 'C1' not in get_batch().job_by_stage
    assert 'C2' not in get_batch().job_by_stage


def test_raises_workflow_error_when_output_from_parent_stage_is_missing(mocker: MockFixture, tmp_path):
    """
    A -> B1 -> C1
    A -> B2 -> C2

    only_stages = [B1]

    Will raise WorkflowError: A is required, but is skipped since outputs for A
    have not been generated.
    """
    conf = create_config(tmp_path, allow_missing_outputs_for_stages=[])
    set_config(conf, tmp_path / 'config.toml')

    with pytest.raises(WorkflowError, match='A: stage is required, but is skipped'):
        run_workflow(mocker, stages=[C1, C2])
