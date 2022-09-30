#!/usr/bin/env python3

"""
Entry point to run the workflow.
"""
import os

import coloredlogs
import click
from cpg_utils.workflows.workflow import get_workflow
from cpg_utils.config import set_config_paths
from stages.multiqc import GvcfMultiQC, CramMultiQC
from stages.seqr_loader import MtToEs

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='INFO', fmt=fmt)


@click.command()
@click.argument('config_paths', nargs=-1)
def main(config_paths: list[str]):
    """
    Run a workflow, using CONFIG_PATHS in the order specified, overriding
    $CPG_CONFIG_PATH if specified.
    """
    if _cpg_config_path_env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths = _cpg_config_path_env_var.split(',') + list(config_paths)
    set_config_paths(list(config_paths))
    get_workflow().run(stages=[MtToEs, CramMultiQC, GvcfMultiQC])


if __name__ == '__main__':
    main()
