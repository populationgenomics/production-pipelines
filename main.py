#!/usr/bin/env python3

"""
Entry point to run the workflow.
"""
import os

import click
import coloredlogs

from cpg_utils import to_path
from cpg_utils.config import set_config_paths, get_config
from cpg_utils.workflows.workflow import run_workflow
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
    base_config_path = to_path(__file__).parent / 'configs' / 'seqr.toml'
    assert base_config_path.exists(), base_config_path
    config_paths = [str(base_config_path)] + list(config_paths)
    if _env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths += _env_var.split(',') + list(config_paths)
    set_config_paths(list(config_paths))

    stages = [GvcfMultiQC, CramMultiQC, MtToEs]
    if last_stages := get_config()['workflow'].get('last_stages'):
        stages = [s for s in stages if s.__name__ in last_stages]
    if only_stages := get_config()['workflow'].get('only_stages'):
        stages = [s for s in stages if s.__name__ in only_stages]
    run_workflow(stages=stages)


if __name__ == '__main__':
    main()
