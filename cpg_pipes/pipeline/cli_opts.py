"""
Common pipeline command line options for Click.
"""
from enum import Enum
from typing import Callable, Type, TypeVar
import click
import click_config_file
import yaml  # type: ignore

from ..providers import (
    StoragePolicyType,
    StatusReporterType,
    InputProviderType,
)
from ..providers.storage import Namespace


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
            help='For the first (not-skipped) stage, if the input for a target does not '
            'exist, just skip this target instead of failing. E.g. if the first '
            'stage is Align, and `sequence.meta` files for a sample do not exist, '
            'remove this sample instead of failing.',
        ),
        click.option(
            '--check-intermediates/--no-check-intermediates',
            'check_intermediates',
            default=True,
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

    # Adding ability to load options from a yaml file
    # using https://pypi.org/project/click-config-file
    def yaml_provider(fp, _):
        """Load options from YAML"""
        with open(fp) as f:
            return yaml.load(f, Loader=yaml.SafeLoader)

    return click_config_file.configuration_option(
        provider=yaml_provider,
        implicit=False,
    )(function)
