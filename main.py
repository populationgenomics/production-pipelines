"""
Entry point.
"""
import logging
import sys
import coloredlogs
import pkg_resources

from cpg_utils.config import get_config
from cpg_utils.workflows.workflow import WorkflowError, Workflow

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='DEBUG', fmt=fmt)


def main():
    if 'stages' not in get_config()['workflow']:
        raise WorkflowError(
            'Please list stages to run with a workflow/stages config parameter'
        )
    
    req_stage_names = get_config()['workflow']['stages']
    if isinstance(req_stage_names, str):
        req_stage_names = [req_stage_names]
    if not req_stage_names:
        raise WorkflowError(
            'Empty list of stages in workflow/stages config parameter, please '
            'specify some stages to run'
        )
    
    avail_stage_by_name = dict()
    for entry_point in pkg_resources.iter_entry_points('stages'):
        class_name = str(entry_point).split("=")[0].strip()
        avail_stage_by_name[class_name] = entry_point
    
    # Check that we were able to find some stages.
    # If not, package probably hasn't been installed properly.
    # Note: Can't use logger here, not yet initiated.
    if len(avail_stage_by_name) == 0:
        logging.error(
            'Error: could not load stages. Has the package been installed?'
            'Please install with pip: (pip install .)',
        )
        sys.exit(1)
    
    # Get the list of stages we want to run, in the order that we want them.
    run_stage_names = [stg for stg in req_stage_names if stg in avail_stage_by_name.keys()]
    if not run_stage_names:
        logging.error(
            f'No known stages requested. Requested: {", ".join(req_stage_names)}, '
            f'available: {", ".join(avail_stage_by_name.keys())}'
        )
    
    logging.debug(f'Running stages: {", ".join(run_stage_names)}')
    stage_classes = [avail_stage_by_name[name].load() for name in run_stage_names]
    workflow = Workflow(stages=stage_classes)
    workflow.run()
    return workflow


if __name__ == '__main__':
    main()
