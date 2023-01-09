"""
Script to remove intermediate file from a workflow output location.
"""

from cpg_workflows import get_workflow
import click


@click.command(no_args_is_help=True)
@click.argument('suffix')
def main(suffix: str):
    """
    Remove arbitrary file from a workflow run output prefix.
    """

    wfl = get_workflow()
    path = wfl.prefix / suffix
    print(f'Using workflow prefix {wfl.prefix}, full path would be {path}')
    if not path.exists():
        print(f'Path {path} does not exist')
    else:
        path.unlink()
        print(f'Removed object {path}')
