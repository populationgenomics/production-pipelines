"""
Create a pipeline from command line options.
"""
import functools
import logging
import os
import tempfile
from typing import Callable, Any, Optional

import click
import toml
import yaml  # type: ignore
from cpg_utils.config import update_dict, write_config, get_config

from .pipeline import Pipeline
from .. import to_path, get_package_path
from ..providers.cpg import analysis_runner_env
from ..utils import exists

logger = logging.getLogger(__file__)

# Fail if a configuration file has unknown parameters.
CLI_CONFIG_FILE_STRICT = True


def file_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    extra_file_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    @param ext: check that the path has the expected extension
    @param must_exist: check that the input file/object/directory exists
    @param extra_file_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. file.vcf.gz and file.vcf.gz.tbi)
    @return: a callback suitable for Click parameter initialization
    """

    def _callback(_: click.Context, param: click.Option, value: Any):
        if not value:
            return value
        vals = value if isinstance(value, tuple) else [value]
        for val in vals:
            if ext:
                val = val.rstrip('/')
                if not val.endswith(f'.{ext}'):
                    raise click.BadParameter(
                        f'Value to {param.name} is expected to have '
                        f'an extension .{ext}, got: {val}'
                    )
            if must_exist:
                if not exists(val):
                    raise click.BadParameter(
                        f"Value to {param.name} {val} doesn't exist or incomplete"
                    )
                if extra_file_suffix:
                    extra_file_path = os.path.splitext(val)[0] + extra_file_suffix
                    if not exists(extra_file_path):
                        raise click.BadParameter(
                            f"An accompanying file {extra_file_path} doesn't "
                            f'exist'
                        )
        return vals
    return _callback


MainFunType = Callable[[Pipeline], None]
WrappedMainFunType = Callable[[], None]


def pipeline_entry_point(
    main_fun: Optional[MainFunType] = None,
    *,
    name: str | None = None,
    description: str | None = None,
):
    """
    Decorator that main function. It creates a Pipeline object and passed it
    to the function. This, input type is a function that takes a Pipeline object: 
    MainFunType; and the output type is a function that takes zero arguments, i.e.
    entry point: WrappedMainFunType.
    
    Also adds an optional `--config` command line option that points to a TOML file,
    to override values in CPG_CONFIG_PATH. For config, see the template with defaults
    config-template.toml in this package.

    >>> @pipeline_entry_point(name='my pipeline')
    >>> def main(pipeline: Pipeline):
    >>>     pipeline.run(stages=[...])
    """
    def _decorator(_main_fun: MainFunType) -> WrappedMainFunType:
        @functools.wraps(_main_fun)
        def _wrappped_main() -> None:
            # Load and puts the default config underneath the user provided config set 
            # by CPG_CONFIG_PATH.
            with to_path(get_package_path() / 'config-template.toml').open() as f:
                config = toml.load(f)
            update_dict(config, get_config())
            tmp_dir = config['workflow'].setdefault('local_tmp_dir', tempfile.mkdtemp())
            write_config(config, tmp_dir)

            # Add analysis-runner environment for the dataset if needed.
            with analysis_runner_env():
                return _main_fun(Pipeline(name, description))

        return _wrappped_main

    if main_fun is None:
        return _decorator
    else:
        return _decorator(main_fun)
