import cpg_utils.config
import pytest

import cpg_workflows.batch
import cpg_workflows.inputs
import cpg_workflows.metamist
import cpg_workflows.stages.gatk_sv
import cpg_workflows.workflow


@pytest.fixture(autouse=True, scope='function')
def pre_and_post_test():
    yield

    # Reset config paths to defaults
    cpg_utils.config.set_config_paths([])

    # Clear pre-existing state before running a new workflow. Must use setattr
    # for this to work so ignore flake8 B010.
    setattr(cpg_utils.config, '_config', None)  # noqa: flake8
    setattr(cpg_workflows.batch, '_batch', None)  # noqa: B010
    setattr(cpg_workflows.workflow, '_workflow', None)  # noqa: B010
    setattr(cpg_workflows.inputs, '_cohort', None)  # noqa: B010
    setattr(cpg_workflows.stages.gatk_sv, '_FASTA', None)  # noqa: B010
    setattr(cpg_workflows.metamist, '_metamist', None)  # noqa: B010
