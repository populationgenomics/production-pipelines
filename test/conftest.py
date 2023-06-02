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

    # Clear pre-existing state before running a new workflow
    setattr(cpg_utils.config, '_config', None)
    setattr(cpg_workflows.batch, '_batch', None)
    setattr(cpg_workflows.workflow, '_workflow', None)
    setattr(cpg_workflows.inputs, '_cohort', None)
    setattr(cpg_workflows.stages.gatk_sv, '_FASTA', None)
    setattr(cpg_workflows.metamist, '_metamist', None)
