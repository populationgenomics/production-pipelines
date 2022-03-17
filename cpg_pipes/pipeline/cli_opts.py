"""
Common pipeline command line options for "click".
"""
from typing import Callable
import click
import click_config_file
import yaml  # type: ignore

from cpg_pipes.storage import Namespace


def pipeline_click_options(function: Callable) -> Callable:
    """
    Decorator to use with click when writing a script that implements a pipeline.
    For example:

    @click.command()
    @click.argument('--custom-argument')
    @pipeline_click_options
    def main(custom_argument: str):
        pass
    """
    options = [
        click.option(
            '-n',
            '--namespace',
            'namespace',
            type=click.Choice([n.lower() for n in Namespace.__members__]),
            callback=lambda c, p, v: getattr(Namespace, v.upper()) if v else None,
            help='The bucket namespace to write the results to',
        ),
        click.option(
            '--analysis-dataset',
            'analysis_dataset',
            help='Dataset name to write the intermediate/joint-calling analysis '
                 'entries in the sample metadata database',
        ),
        click.option(
            '--input-dataset',
            'input_datasets',
            multiple=True,
            help='Only read samples that belong to the dataset(s). '
                 'Can be set multiple times.',
        ),
        click.option(
            '--ped-file',
            'ped_files',
            multiple=True,
            help='PED file (will override sample-meatadata family data if available)'
        ),
        click.option(
            '--first-stage',
            'first_stage',
            help='Skip previous stages and pick their expected results if further '
                 'stages depend on thems',
        ),
        click.option(
            '--last-stage',
            'last_stage',
            help='Finish the pipeline after this stage',
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
            help='Only take these samples (can be set multiple times)',
        ),
        click.option(
            '--force-sample',
            'force_samples',
            multiple=True,
            help='Force reprocessing these samples. Can be set multiple times.',
        ),
        click.option(
            '--version', '--output-version'
            'version',
            type=str,
            help='Pipeline version. Default is a timestamp',
        ),
        click.option(
            '--keep-scratch/--remove-scratch', 
            'keep_scratch', 
            default=False,
            is_flag=True,
        ),
        click.option('--dry-run', 'dry_run', is_flag=True),
        click.option(
            '--check-smdb-seq-existence/--no-check-smdb-seq-existence',
            'check_smdb_seq',
            default=False,
            is_flag=True,
            help='Check that files in sequence.meta exist'
        ),
        click.option(
            '--skip-samples-with-missing-input',
            'skip_samples_with_missing_input',
            default=False,
            is_flag=True,
            help='For the first (not-skipped) stage, if the input for a target does not'
                 'exist, just skip this target instead of failing. E.g. if the first'
                 'stage is CramStage, and sequence.meta files for a sample do not exist,'
                 'remove this sample instead of failing.'
        ),
        click.option(
            '--check-intermediate-existence/--no-check-intermediate-existence',
            'check_intermediates',
            default=True,
            is_flag=True,
            help='Within jobs, check all in-job intermediate files for possible reuse. '
                 'If set to False, will overwrite all intermediates. '
        ),
        click.option(
            '--check-job-expected-outputs-existence/--no-check-job-expected-outputs-existence',
            'check_expected_outputs',
            default=True,
            is_flag=True,
            help='Before running a job, check if its input already exists. '
                 'If it exists, submit a [reuse] job instead. '
                 'Works nicely with --previous-batch-tsv/--previous-batch-id options.',
        ),
        click.option(
            '--update-smdb-analyses/--no-update-smdb-analyses',
            'update_smdb_analyses',
            is_flag=True,
            default=True,
            help='Create analysis entries for queued/running/completed jobs'
        ),
        click.option(
            '--validate-smdb-analyses/--no-validate-smdb-analyses',
            'validate_smdb_analyses',
            is_flag=True,
            default=False,
            help='Validate existing analysis entries by checking if a.output exists on '
                 'the bucket. Set the analysis entry to "failure" if output doesn\'t '
                 'exist'
        ),
        click.option(
            '--previous-batch-tsv',
            'previous_batch_tsv_path',
            help='A list of previous successful attempts from another batch, dumped '
                 'from from the Batch database (the "jobs" table joined on '
                 '"job_attributes"). If the intermediate output for a job exists in '
                 'a previous attempt, it will be passed forward, and a [reuse] job will '
                 'be submitted.'
        ),
        click.option(
            '--previous-batch-id',
            'previous_batch_id',
            help='6-letter ID of the previous successful batch (corresponds to the '
                 'directory name in the batch logs. e.g. feb0e9 in '
                 'gs://cpg-seqr-main-tmp/hail/batch/feb0e9'
        ),
        click.option(
            '--local-dir',
            'local_dir',
        ),
    ]
    # Click shows options in a reverse order, so inverting the list back:
    for opt in options[::-1]:
        function = opt(function)
    
    # Add ability to load options from a yaml file
    # using https://pypi.org/project/click-config-file/
    def yaml_provider(fp, _):
        """Load options from YAML"""
        with open(fp) as f:
            return yaml.load(f, Loader=yaml.SafeLoader)
    function = click_config_file.configuration_option(
        provider=yaml_provider
    )(function)
 
    return function
