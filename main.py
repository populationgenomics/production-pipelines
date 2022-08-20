"""
Entry point.
"""
import logging
import sys
import coloredlogs
import pkg_resources

from cpg_utils.config import get_config
from cpg_utils.flows.workflow import WorkflowError, Workflow

logger = logging.getLogger(__file__)

coloredlogs.install(
    level='DEBUG', fmt='%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
)


if 'stages' not in get_config()['workflow']:
    raise WorkflowError(
        'Please list stages to run with a workflow/stages config parameter'
    )

requested_stages = get_config()['workflow']['stages']
if not requested_stages:
    raise WorkflowError(
        'Empty list of stages in workflow/stages config parameter, please '
        'specify some stages to run'
    )

avail_stages = dict()
for entry_point in pkg_resources.iter_entry_points('stages'):
    class_name = str(entry_point).split("=")[0].strip()
    avail_stages[class_name] = entry_point

# Check that we were able to find some stages.
# If not, package probably hasn't been installed properly.
# Note: Can't use logger here, not yet initiated.
if len(avail_stages) == 0:
    logger.error(
        'Error: could not load stages. Has the package been installed?'
        'Please install with pip: (pip install .)',
    )
    sys.exit(1)

# Get the list of stages we want to run, in the order that we want them.
run_stages = [_stage for _stage in requested_stages if _stage in avail_stages.keys()]
if not run_stages:
    logger.error(
        f'No known stages requested. Requested: {", ".join(requested_stages)}, '
        f'availabel: {", ".join(avail_stages.keys())}'
    )

logger.debug(f'Running stages: {", ".join(run_stages)}')

workflow = Workflow(stages=[stage.load() for stage in avail_stages.values()])
workflow.run()
