#!/usr/bin/env python3

"""
Entry point to run the workflow.
"""
import os

import click
import coloredlogs

from cpg_utils import to_path
from cpg_utils.config import set_config_paths
from cpg_utils.workflows.workflow import run_workflow
from cpg_workflows.stages.large_cohort import LoadVqsr, Frequencies
from cpg_workflows.stages.multiqc import GvcfMultiQC, CramMultiQC
from cpg_workflows.stages.fastqc import FastQCMultiQC
from cpg_workflows.stages.seqr_loader import MtToEs

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='INFO', fmt=fmt)


WORKFLOWS = {
    'pre_alignment': [FastQCMultiQC],
    'seqr_loader': [MtToEs, GvcfMultiQC, CramMultiQC],
    'large_cohort': [LoadVqsr, Frequencies, GvcfMultiQC, CramMultiQC],
}


@click.command()
@click.argument('workflow', type=click.Choice(list(WORKFLOWS.keys())))
@click.option('--config', 'config_paths', multiple=True)
def main(workflow: str, config_paths: list[str]):
    """
    Run a workflow, using CONFIG_PATHS in the order specified, overriding
    $CPG_CONFIG_PATH if specified.
    """
    base_config_path = (
        to_path(__file__).parent / 'configs' / 'defaults' / f'{workflow}.toml'
    )
    assert base_config_path.exists(), base_config_path
    config_paths = [str(base_config_path)] + list(config_paths)
    if _env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths += _env_var.split(',') + list(config_paths)
    set_config_paths(list(config_paths))

    run_workflow(stages=WORKFLOWS[workflow])


if __name__ == '__main__':
    main()
