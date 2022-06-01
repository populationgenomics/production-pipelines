"""
Create a pipeline from command line options.
"""

import logging
import os
import tempfile
import uuid
from enum import Enum
from typing import Callable, Type, TypeVar, Any

import click
import toml
import yaml  # type: ignore
from click import option
from cpg_utils.config import set_config_path
from cpg_utils.hail_batch import Namespace, Path

from .exceptions import PipelineError
from .pipeline import Pipeline
from .. import to_path
from ..providers import (
    StatusReporterType,
    InputProviderType,
)
from ..providers.cpg import build_infra_config, set_config
from ..providers.cpg.images import CpgImages
from ..providers.cpg.inputs import SmdbInputProvider
from ..providers.cpg.refdata import CpgRefData
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


def file_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    accompanying_metadata_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    @param ext: check that the path has the expected extension
    @param must_exist: check that the input file/object/directory exists
    @param accompanying_metadata_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. genomes.mt and genomes.metadata.ht)
    @return: a callback suitable for Click parameter initialization
    """

    def _callback(_: click.Context, param: click.Option, value: Any):
        if value is None:
            return value
        if ext:
            assert isinstance(value, str), value
            value = value.rstrip('/')
            if not value.endswith(f'.{ext}'):
                raise click.BadParameter(
                    f'The argument {param.name} is expected to have '
                    f'an extension .{ext}, got: {value}'
                )
        if must_exist:
            if not exists(value):
                raise click.BadParameter(f"{value} doesn't exist or incomplete")
            if accompanying_metadata_suffix:
                accompanying_metadata_fpath = (
                    os.path.splitext(value)[0] + accompanying_metadata_suffix
                )
                if not exists(accompanying_metadata_fpath):
                    raise click.BadParameter(
                        f"An accompanying file {accompanying_metadata_fpath} doesn't "
                        f'exist'
                    )
        return value
    return _callback


class Options:
    """List of Click options, to allow access from attributes."""
    namespace: option(
        '-n',
        '--namespace',
        type=choice_from_enum(Namespace),
        callback=val_to_enum(Namespace),
        required=True,
        help='The bucket namespace to write the results to',
    )
    analysis_dataset: option(
        '--analysis-dataset',
        required=True,
        help='Dataset name to write cohort and pipeline level intermediate files',
    )
    datasets: option(
        '--dataset', 
        'datasets',
        multiple=True,
        help='Only read samples that belong to the given dataset(s). '
        'Can be set multiple times.',
    )
    first_stage: option(
        '--first-stage',
        help='Skip previous stages and pick their expected results if further '
        'stages depend on them',
    )
    last_stage: option(
        '--last-stage',
        help='Finish the pipeline after this stage',
    )
    skip_datasets: option(
        '--skip-dataset',
        'skip_datasets',
        multiple=True,
        help='Don\'t process specified datasets. Can be set multiple times.',
    )
    skip_samples: option(
        '--skip-sample',
        '-S',
        'skip_samples',
        multiple=True,
        help='Don\'t process specified samples. Can be set multiple times.',
    )
    only_samples: option(
        '--only-sample',
        '-s',
        'only_samples',
        multiple=True,
        help='Only process these samples (can be set multiple times)',
    )
    force_samples: option(
        '--force-sample',
        'force_samples',
        multiple=True,
        help='Force reprocessing these samples. Can be set multiple times.',
    )
    sequencing_type: option(
        '--sequencing-type',
        type=choice_from_enum(SequencingType),
        callback=val_to_enum(SequencingType),
        help='Limit to data with this sequencing type',
        default=SequencingType.GENOME.value,
    )
    name: option(
        '--name',
        help='Name of the pipeline (to prefix output paths)'
    )
    description: option(
        '--description',
        help='Description of the pipeline (to appear in the Batch GUI)'
    )
    version: option(
        '--version',
        '--output-version',
        'version',
        type=str,
        help='Pipeline version. Default is a timestamp',
    )
    status_reporter: option(
        '--status-reporter',
        type=choice_from_enum(StatusReporterType),
        callback=val_to_enum(StatusReporterType),
        help='Report jobs statuses and results',
    )
    input_provider: option(
        '--input-provider',
        type=choice_from_enum(InputProviderType),
        callback=val_to_enum(InputProviderType),
        default=InputProviderType.CPG.value,
        help=f'Source of inputs. '
        f'For "--input-source={InputProviderType.CSV.value}", '
        f'use --input-csv to specify a CSV file location',
    )
    input_csv: option(
        '--input-csv',
        callback=file_validation_callback(ext='csv', must_exist=True),
        help=f'CSV file to provide with --input-provider={InputProviderType.CSV.value}',
    )
    keep_scratch: option(
        '--keep-scratch/--remove-scratch',
        'keep_scratch',
        default=False,
        is_flag=True,
    )
    dry_run: option(
        '--dry-run', 
        is_flag=True, 
        help='Do not actually submit Batch, but only print jobs commands to '
             'stdout. Essencially just tell to call `batch.run(..., dry_run=True)`'
    )
    skip_samples_with_missing_input: option(
        '--skip-samples-with-missing-input',
        default=False,
        is_flag=True,
        help='For the first (not-skipped) stage, if the input for a target does '
        'not exist, just skip this target instead of failing. E.g. if the first '
        'stage is Align, and `sample.alignment_input` for a sample do not exist, '
        'remove this sample, instead of failing. In order words, ignore samples '
        'that are missing results from skipped stages.',
    )
    check_intermediates: option(
        '--check-intermediates/--no-check-intermediates',
        'check_intermediates',
        default=False,
        is_flag=True,
        help='Within jobs, check all in-job intermediate files for possible reuse. '
        'If set to False, will overwrite all intermediates.',
    )
    check_inputs: option(
        '--check-inputs/--no-check-inputs',
        'check_inputs',
        default=False,
        is_flag=True,
        help='Check input file existence (e.g. FASTQ files). If they are missing '
        'the --skip-samples-with-missing-input option controls whether such '
        'should be ignored, or raise an error.',
    )
    check_expected_outputs: option(
        '--check-expected-outputs/--no-check-expected-outputs',
        'check_expected_outputs',
        default=True,
        is_flag=True,
        help='Before running a stage, check if its input already exists. '
        'If it exists, submit a [reuse] job instead. '
        'Works nicely with --previous-batch-tsv/--previous-batch-id options.',
    )
    hail_pool_label: option(
        '--hail-pool-label',
        help='Private pool label. Would submit batches to a Hail Batch private '
             'pool with this label'
    )
    local_tmp_dir: option(
        '--local-dir',
        'local_tmp_dir',
        default=tempfile.mkdtemp(),
        type=click.Path(exists=True, dir_okay=True, file_okay=False),
        help='Local directory for temporary files. Usually takes a few kB. '
        'If not provided, a temp folder will be created',
    )
    slack_channel: option(
        '--slack-channel',
        help='Slack channel to send status reports with CPG status reporter.',
    )
    smdb_errors_are_fatal: option(
        '--smdb-errors-are-fatal/--no-smdb-errors-are-fatal',
        'smdb_errors_are_fatal',
        default=True,
        help='When the sample-metadata database API returns an error from the '
             'database, only show the error and continue, instead of crashing',
    )
    skip_samples_stages: option(
        '--skip-samples-stages',
        type=dict,
        help='Map of stages to lists of samples, to skip for specific stages.',
    )
    image_registry_prefix: option(
        '--image-registry-prefix',
        help='Prefix to find docker images referenced in `image_config_yaml_path`',
    )
    image_config: option(
        '--image-config',
        'image_config_yaml_path',
        help='Path to a YAML file containing - for each image name used in the '
             'pipelines - a URL of a corresponding docker image, optionally '
             'relative to the prefix specified by `image_registry_prefix`.',
    )
    reference_prefix: option(
        '--reference-prefix',
        help='Prefix to find reference files',
    )
    gcp_project: option(
        '--gcp-project',
        help='GCP project',
    )
    hail_billing_project: option(
        '--hail-billing-project',
        help='Hail billing project',
    )
    hail_bucket: option(
        '--hail-bucket',
        help='Hail bucket',
    )
    web_url_template: option(
        '--web-url-template',
        help='Template to build HTTP URLs matching the dataset_path of category '
             '"web". Should be parametrised by namespace and dataset in Jinja format',
    )

    def __init__(self, kwargs: dict):
        self.remaining_kwargs = kwargs.copy()
        for k, v in kwargs.items():
            if k in self.__annotations__:
                setattr(self, k, v)
                del self.remaining_kwargs[k]


def create_pipeline(**kwargs) -> 'Pipeline':
    """
    Create a Pipeline instance. All options correspond to command line parameters
    described in `pipeline_click_options` in the `cli_opts` module
    """
    options = Options(kwargs)

    if not os.getenv('CPG_CONFIG_PATH'):
        infra_config = build_infra_config(
            dataset=options.analysis_dataset, 
            namespace=options.namespace,
            gcp_project=options.gcp_project,
            image_registry_prefix=options.image_registry_prefix,
            reference_prefix=options.reference_prefix,
            web_url_template=options.web_url_template,
        )
        set_config(infra_config, to_path(options.local_tmp_dir))

    refs = CpgRefData()
    images = CpgImages()
    status_reporter: StatusReporter | None = None
    input_provider: InputProvider | None = None
    
    analysis_dataset = Dataset(
        options.analysis_dataset, 
        namespace=options.namespace
    )

    if (
        options.input_provider == InputProviderType.CPG or 
        options.status_reporter == StatusReporterType.CPG
    ):
        smdb = SMDB(analysis_dataset.name)
        if options.status_reporter == StatusReporterType.CPG:
            status_reporter = CpgStatusReporter(
                smdb=smdb,
                images=images,
                slack_channel=options.slack_channel,
            )
        if options.input_provider == InputProviderType.CPG:
            input_provider = SmdbInputProvider(
                smdb,
                smdb_errors_are_fatal=options.smdb_errors_are_fatal,
            )

    if options.input_provider == InputProviderType.CSV:
        if not options.input_csv:
            raise PipelineError(
                f'input_csv (--input-csv) should be provided '
                f'with input_provider_type=InputProviderType.CSV '
                f'(--input-provider {InputProviderType.CSV.value})'
            )
        input_provider = CsvInputProvider(to_path(options.input_csv).open())

    return Pipeline(
        namespace=options.namespace,
        name=options.name or analysis_dataset.name,
        description=options.description,
        analysis_dataset_name=options.analysis_dataset,
        refs=refs,
        images=images,
        input_provider=input_provider,
        sequencing_type=options.sequencing_type,
        status_reporter=status_reporter,
        datasets=options.datasets,
        skip_datasets=options.skip_datasets,
        skip_samples=options.skip_samples,
        only_samples=options.only_samples,
        force_samples=options.force_samples,
        first_stage=options.first_stage,
        last_stage=options.last_stage,
        version=options.version,
        check_inputs=options.check_inputs,
        check_intermediates=options.check_intermediates,
        check_expected_outputs=options.check_expected_outputs,
        skip_samples_with_missing_input=options.skip_samples_with_missing_input,
        local_dir=options.local_tmp_dir,
        dry_run=options.dry_run,
        keep_scratch=options.keep_scratch,
        skip_samples_stages=options.skip_samples_stages,
        hail_pool_label=options.hail_pool_label,
        **options.remaining_kwargs,
    )


def pipeline_options(function: Callable) -> Callable:
    """
    Decorator to use with Click, that adds common pipeline options.
    Useful to use in a script that creates a pipeline. Arguments can
    be passed directly to `create_pipeline`, for example:

    >>> @click.command()
    >>> @pipeline_options
    >>> def main(**kwargs):
    >>>     pipeline = create_pipeline(**kwargs)

    Domain specific options can be added using more click decorators before
    `@pipeline_options`, e.g.:

    >>> @click.command()
    >>> @click.argument('--custom-argument')
    >>> @pipeline_options
    >>> def main(custom_argument: str, **kwargs):
    >>>     pipeline = create_pipeline(**kwargs, config=dict(myarg=custom_argument))
    """

    # Applying decorators. Doing that in reverse order, because Click actually
    # inverts the order of shown options, assuming the decorators order of
    # application which is bottom to top.
    for opt in list(Options.__annotations__.values())[::-1]:
        function = opt(function)

    # Adding a special option to read defaults from a user-provided configuration file. 
    # Adding it last, so we are able to compare parameters in the file with the "normal" 
    # options that would be populated in function.__click_params__ by click at this 
    # point.
    return option(
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
        if not yaml_path:
            return
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
            if CLI_CONFIG_FILE_STRICT:
                raise click.BadOptionUsage(param, msg, ctx)
            else:
                logger.warning(msg)

        ctx.default_map = ctx.default_map or {}
        ctx.default_map.update(d)
    return config_callback
