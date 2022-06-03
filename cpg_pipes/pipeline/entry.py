"""
Create a pipeline from command line options.
"""
import functools
import logging
import os
import sys
import tempfile
from enum import Enum
from typing import Callable, Type, TypeVar, Any, Optional, Union

import click
import toml
import yaml  # type: ignore
from cpg_utils.config import get_config, update_dict, set_config_path
from cpg_utils.hail_batch import Namespace
from slugify import slugify

from .exceptions import PipelineError
from .pipeline import Pipeline, StageDecorator
from .. import to_path, get_package_path
from ..jobs import haplotype_caller, joint_genotyping, vep, align
from ..providers import (
    StatusReporterType,
    InputProviderType,
)
from ..providers.cpg import complete_infra_config, overwrite_config
from ..providers.cpg.inputs import CpgInputProvider
from ..providers.cpg.smdb import SMDB
from ..providers.cpg.status import CpgStatusReporter
from ..providers.inputs import InputProvider, CsvInputProvider
from ..providers.status import StatusReporter
from ..targets import Dataset
from ..types import SequencingType
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
        @click.command()
        @click.option(
            '--config',
            'config_paths',
            multiple=True,
            callback=file_validation_callback(ext='toml', must_exist=True),
            help='Read configuration from a TOML files. Can be set multiple times.',
        )
        @functools.wraps(_main_fun)
        def _wrappped_main(config_paths: list[str]) -> None:
            with to_path(get_package_path() / 'config-template.toml').open() as f:
                config = toml.load(f)
            if os.environ.get('CPG_CONFIG_PATH'):
                update_dict(config, get_config())
            for path in config_paths:
                with to_path(path).open() as f:
                    update_dict(config, toml.load(f)) 
            config = complete_infra_config(config)
            # Add contents of CPG_CONFIG_PATH on top, write the result 
            # and reset CPG_CONFIG_PATH:
            overwrite_config(config)
            pipeline = Pipeline(name=name, description=description)
            return _main_fun(pipeline)

        return _wrappped_main

    if main_fun is None:
        return _decorator
    else:
        return _decorator(main_fun)
