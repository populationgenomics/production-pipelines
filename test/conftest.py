import os
from unittest import mock

import pytest
from google.auth import environment_vars

import cpg_utils.config
import cpg_workflows.batch
import cpg_workflows.inputs
import cpg_workflows.metamist
import cpg_workflows.stages.gatk_sv
import cpg_workflows.workflow


@pytest.fixture(autouse=True, scope='function')
def pre_and_post_test():
    # Set a dummy google cloud project to avoid errors when running tests for tests
    # that use the google cloud.
    with mock.patch.dict(
        os.environ,
        {environment_vars.PROJECT: 'dummy-project-for-tests'},
    ):
        yield

    # Reset config paths to defaults
    cpg_utils.config.set_config_paths([])

    # Clear pre-existing state before running a new workflow. Must use setattr
    # for this to work so ignore flake8 B010.
    setattr(cpg_utils.config, '_config', None)  # noqa: B010
    setattr(cpg_workflows.batch, '_batch', None)  # noqa: B010
    setattr(cpg_workflows.workflow, '_workflow', None)  # noqa: B010
    setattr(cpg_workflows.inputs, '_cohort', None)  # noqa: B010
    setattr(cpg_workflows.stages.gatk_sv, '_FASTA', None)  # noqa: B010
    setattr(cpg_workflows.metamist, '_metamist', None)  # noqa: B010
    setattr(cpg_workflows.inputs, '_multicohort', None)  # noqa: B010
