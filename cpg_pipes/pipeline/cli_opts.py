"""
Common pipeline command line options for Click.
"""
import logging
from enum import Enum
from typing import Callable, Type, TypeVar
import click
import yaml  # type: ignore

from ..providers import (
    StoragePolicyType,
    StatusReporterType,
    InputProviderType,
)
from ..providers.storage import Namespace

logger = logging.getLogger(__file__)

# Whether to fail if a configuration file has unknown parameters.
CONFIG_FILE_STRICT = True


def choice_from_enum(cls: Type[Enum]) -> click.Choice:
    """
    Create click.Choice from Enum items.
    """
    return click.Choice([n.lower() for n in cls.__members__])


T = TypeVar('T', bound=Enum)


def val_to_enum(cls: Type[T]) -> Callable:
    """
    Callback to parse a value into an Enum item.
    """

    def _callback(ctx, param, val: str) -> T | None:
        if val is None:
            return None
        d = {name.lower(): item for name, item in cls.__members__.items()}
        if val not in d:
            raise click.BadParameter(
                f'Available options: {[n.lower() for n in cls.__members__]}'
            )
        return d[val]

    return _callback


def pipeline_click_options(function: Callable) -> Callable:
    """
    Decorator to use with Click that adds common pipeline options.
    Useful to use with a script that implements a pipeline. Arguments can
    be passed directly to `create_pipeline`, for example:

        @click.command()
        @pipeline_click_options
        def main(**kwargs):
            pipeline = create_pipeline(**kwargs)

    New options can be added by adding more click decorators before
    `@pipeline_click_options`, e.g.:

        @click.command()
        @click.argument('--custom-argument')
        @pipeline_click_options
        def main(custom_argument: str):
            pipeline = create_pipeline(**kwargs, config=dict(myarg=custom_argument))
    """

    options = [
        click.option(
            '-n',
            '--namespace',
            'namespace',
            type=choice_from_enum(Namespace),
            callback=val_to_enum(Namespace),
            required=True,
            help='The bucket namespace to write the results to',
        ),
        click.option(
            '--analysis-dataset',
            'analysis_dataset',
            help='Dataset name to write cohort and pipeline level intermediate files',
            required=True,
        ),
        click.option(
            '--dataset',
            'datasets',
            multiple=True,
            help='Only read samples that belong to the given dataset(s). '
            'Can be set multiple times.',
        ),
        click.option(
            '--first-stage',
            'first_stage',
            help='Skip previous stages and pick their expected results if further '
            'stages depend on them',
        ),
        click.option(
            '--last-stage',
            'last_stage',
            help='Finish the pipeline after this stage',
        ),
        click.option(
            '--skip-dataset',
            '-D',
            'skip_datasets',
            multiple=True,
            help='Don\'t process specified datasets. Can be set multiple times.',
        ),
        click.option(
            '--skip-sample',
            '-S',
            'skip_samples',
            multiple=True,
            help='Don\'t process specified samples. Can be set multiple times.',
        ),
        click.option(
            '--only-sample',
            '-s',
            'only_samples',
            multiple=True,
            help='Only process these samples (can be set multiple times)',
        ),
        click.option(
            '--force-sample',
            'force_samples',
            multiple=True,
            help='Force reprocessing these samples. Can be set multiple times.',
        ),
        click.option(
            '--version',
            '--output-version',
            'version',
            type=str,
            help='Pipeline version. Default is a timestamp',
        ),
        click.option(
            '--storage-policy',
            'storage_policy_type',
            type=choice_from_enum(StoragePolicyType),
            callback=val_to_enum(StoragePolicyType),
            default=StoragePolicyType.CPG.value,
            help='Storage policy is used to determine bucket names for intermediate '
            'and output files',
        ),
        click.option(
            '--status-reporter',
            'status_reporter_type',
            type=choice_from_enum(StatusReporterType),
            callback=val_to_enum(StatusReporterType),
            help='Report jobs statuses and results',
        ),
        click.option(
            '--input-provider',
            'input_provider_type',
            type=choice_from_enum(InputProviderType),
            callback=val_to_enum(InputProviderType),
            default=InputProviderType.SMDB.value,
            help=f'Source of inputs. '
            f'For "--input-source={InputProviderType.CSV.value}", '
            f'use --input-csv to specify a CSV file location',
        ),
        click.option(
            '--input-csv',
            'input_csv',
            help=f'CSV file to provide with --input-provider={InputProviderType.CSV.value}',
        ),
        click.option(
            '--keep-scratch/--remove-scratch',
            'keep_scratch',
            default=False,
            is_flag=True,
        ),
        click.option('--dry-run', 'dry_run', is_flag=True),
        click.option(
            '--skip-samples-with-missing-input',
            'skip_samples_with_missing_input',
            default=False,
            is_flag=True,
            help='For the first (not-skipped) stage, if the input for a target does '
            'not exist, just skip this target instead of failing. E.g. if the first '
            'stage is Align, and `sample.alignment_input` for a sample do not exist, '
            'remove this sample, instead of failing. In order words, ignore samples '
            'that are missing results from skipped stages.',
        ),
        click.option(
            '--check-intermediates/--no-check-intermediates',
            'check_intermediates',
            default=False,
            is_flag=True,
            help='Within jobs, check all in-job intermediate files for possible reuse. '
            'If set to False, will overwrite all intermediates.',
        ),
        click.option(
            '--check-inputs/--no-check-inputs',
            'check_inputs',
            default=False,
            is_flag=True,
            help='Check input file existence (e.g. FASTQ files). If they are missing '
            'the --skip-samples-with-missing-input option controls whether such '
            'should be ignored, or raise an error.',
        ),
        click.option(
            '--check-expected-outputs/--no-check-expected-outputs',
            'check_expected_outputs',
            default=True,
            is_flag=True,
            help='Before running a stage, check if its input already exists. '
            'If it exists, submit a [reuse] job instead. '
            'Works nicely with --previous-batch-tsv/--previous-batch-id options.',
        ),
        click.option(
            '--local-dir',
            'local_dir',
            help='Local directory for temporary files. Usually takes a few kB. '
            'If not provided, a temp folder will be created',
        ),
        click.option(
            '--slack-channel',
            'slack_channel',
            help='Slack channel to send status reports with CPG status reporter.',
        ),
        click.option(
            '--smdb-errors-are-fatal/--no-smdb-errors-are-fatal',
            'smdb_errors_are_fatal',
            default=True,
            help='When the sample-metadata database API returns an error from the '
                 'database, only show the error and continue, instead of crashing',
        ),
        click.option(
            '--skip-samples-stages',
            'skip_samples_stages',
            type=dict,
            help='Map of stages to lists of samples, to skip for specific stages.',
        ),
    ]

    # Applying decorators. Doing that in reverse order, because Click actually
    # inverts the order of shown options, assuming the decorators order of
    # application which is bottom to top.
    for opt in options[::-1]:
        function = opt(function)
    
    # Adding a special option to read defaults from a user-provided configuration file. 
    # Adding it last, so we are able to compare parameters in the file with the "normal" 
    # options that would be populated in function.__click_params__ by click at this 
    # point.
    return click.option(
        '--config',
        'config',
        type=click.Path(exists=True, dir_okay=False, readable=True),
        help='Read configuration from a YAML FILE.',
        callback=get_config_callback(getattr(function, '__click_params__', [])),
        expose_value=False,  # don't poss parameter to `main()`
        is_eager=True,  # process before other options, to make sure we don't fail 
                        # is required command line options are missing
    )(function)


def get_config_callback(defined_options: list[click.Option]):
    """
    Create callback for the --config option. We want to be able to check what params 
    are defined by "normal" click options to compare them with the params in the 
    config file, so for this reason we parametrise the callback with `defined_options`
    by wrapping it with another function `get_config_callback`.
    """
    def config_callback(ctx, param, yaml_path):
        """Read default option values from a YAML config file"""
        try:
            with open(yaml_path) as f:
                d = yaml.load(f, Loader=yaml.SafeLoader)
        except Exception as e:
            raise click.BadOptionUsage(
                param, f'Error reading configuration file: {e}', ctx
            )
        
        defined_opts = [_opt.name for _opt in defined_options]
        unused_defaults = [k for k in d.keys() if k not in defined_opts]
        if unused_defaults:
            msg = (
                f'Found unknown option(s) in the config {yaml_path}: '
                f'{unused_defaults}'
            )
            if CONFIG_FILE_STRICT:
                raise click.BadOptionUsage(param, msg, ctx)
            else:
                logger.warning(msg)

        ctx.default_map = ctx.default_map or {}
        ctx.default_map.update(d)
    return config_callback
