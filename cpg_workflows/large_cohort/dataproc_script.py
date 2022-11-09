#!/usr/bin/env python

"""
Generic script to run a function on dataproc.
"""

import click
from cpg_utils.hail_batch import start_query_context
from cpg_utils import to_path
from importlib import import_module


@click.command()
@click.argument('import_module_name')
@click.argument('function_name')
@click.argument('function_str_args', nargs=-1)
@click.option('-p', '--path', 'function_path_args', multiple=True)
def main(
    import_module_name: str,
    function_name: str,
    function_str_args: list[str],
    function_path_args: list[str],
):
    module = import_module(import_module_name)
    func = getattr(module, function_name)
    start_query_context()
    func(*[to_path(path) for path in function_path_args], *function_str_args)


if __name__ == '__main__':
    main()
