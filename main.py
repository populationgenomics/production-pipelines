#!/usr/bin/env python3

"""
Entry point to run workflows.
"""
import os
import click
import coloredlogs

from cpg_utils import to_path
from cpg_utils.config import set_config_paths
from cpg_workflows import defaults_config_path
from cpg_workflows.workflow import run_workflow

# Importing WORKFLOWS from the correct file based on the environment:
SM_ENVIRONMENT = os.getenv('PL_ENVIRONMENT', 'production').lower()
if SM_ENVIRONMENT == 'production':
    from workflows_prd import WORKFLOWS
else:
    from workflows_dev import WORKFLOWS


@click.command(no_args_is_help=True)
@click.argument('workflow', required=False)
@click.option(
    '--config',
    'config_paths',
    multiple=True,
    help='Add configuration files to the files specified $CPG_CONFIG_PATH.'
    'Configs are merged left to right, meaning the rightmost file has the'
    'highest priority.',
)
@click.option(
    '--list-workflows',
    'list_workflows',
    is_flag=True,
    help='Only list possible values for WORKFLOW (and available last stages)',
)
@click.option(
    '--list-last-stages',
    'list_last_stages',
    is_flag=True,
    help='Only list possible end stages for a workflow, that can be specified with `workflow/last_stages` in config',
)
@click.option(
    '--dry-run',
    'dry_run',
    is_flag=True,
    help='Dry run: do not actually communicate with Metamist or Hail Batch, '
    'instead only print a final config and stages to be run',
)
@click.option(
    '--verbose',
    'verbose',
    is_flag=True,
)
def main(
    workflow: str,
    config_paths: list[str],
    list_workflows: bool,
    list_last_stages: bool,
    dry_run: bool,
    verbose: bool,
):
    """
    Run a Hail Batch workflow specified as a positional command line argument [WORKFLOW]
    """
    fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
    coloredlogs.install(level='DEBUG' if verbose else 'INFO', fmt=fmt)

    if not workflow and not list_workflows:
        click.echo('You must specify WORKFLOW as a first positional command line argument.')
    if not workflow or list_workflows or workflow == 'list':
        click.echo('Available values for WORKFLOW (and corresponding last stages):')
        for wfl, last_stages in WORKFLOWS.items():
            click.echo(f'\t{wfl} ({", ".join(s.__name__ for s in last_stages)})')
        return

    if list_last_stages:
        click.echo(
            f'Available last stages that can be listed with '
            f'workflow/last_stages for the current workflow "{workflow}":',
        )
        click.echo(f'{", ".join(s.__name__ for s in WORKFLOWS[workflow])}')
        return

    wfl_conf_path = to_path(__file__).parent / f'configs/defaults/{workflow}.toml'
    assert wfl_conf_path.exists(), wfl_conf_path

    for path in config_paths:
        assert to_path(path).exists(), path

    config_paths = os.environ['CPG_CONFIG_PATH'].split(',') + list(config_paths)
    # Assuming the defaults is already loaded in __init__.py:
    assert to_path(config_paths[0]) == defaults_config_path
    # Inserting after the "defaults" config, but before user configs:
    set_config_paths(config_paths[:1] + [str(wfl_conf_path)] + config_paths[1:])

    run_workflow(stages=WORKFLOWS[workflow], dry_run=dry_run)


if __name__ == '__main__':
    main()
