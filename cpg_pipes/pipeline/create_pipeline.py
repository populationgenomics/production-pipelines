"""
Provides a `create_pipeline` function that create an instance of a `Pipeline`
using concrete implementations of the input and storage providers.
"""
import logging
from typing import Callable
from hailtop.batch import Batch

from .pipeline import Pipeline, Stage
from .exceptions import PipelineError
from ..targets import Cohort
from .. import Path, to_path
from ..hb.batch import setup_batch, get_billing_project
from ..providers import (
    Cloud,
    Namespace,
    StoragePolicyType,
    StatusReporterType,
    InputProviderType,
)
from ..providers.cpg import (
    CpgStorageProvider,
    SmdbStatusReporter,
    SmdbInputProvider,
)
from ..providers.cpg.smdb import SMDB
from ..providers.inputs import InputProvider, CsvInputProvider
from ..providers.status import StatusReporter
from ..refdata import RefData

logger = logging.getLogger(__file__)


def create_pipeline(
    analysis_dataset: str,
    name: str,
    namespace: Namespace,
    description: str | None = None,
    storage_policy_type: StoragePolicyType = StoragePolicyType.CPG,
    cloud: Cloud = Cloud.GS,
    status_reporter_type: StatusReporterType = None,
    input_provider_type: InputProviderType = InputProviderType.SMDB,
    input_csv: str | None = None,
    stages: list[Callable[..., Stage]] | None = None,
    dry_run: bool = False,
    keep_scratch: bool = True,
    version: str | None = None,
    skip_samples_with_missing_input: bool = False,
    check_intermediates: bool = True,
    check_expected_outputs: bool = True,
    first_stage: str | None = None,
    last_stage: str | None = None,
    config: dict | None = None,
    datasets: list[str] | None = None,
    skip_datasets: list[str] | None = None,
    skip_samples: list[str] | None = None,
    only_samples: list[str] | None = None,
    force_samples: list[str] | None = None,
    ped_files: list[Path] | None = None,
    local_dir: Path | None = None,
) -> 'Pipeline':
    """
    Create a Pipeline instance. All options correspond to command line parameters
    described in `pipeline_click_options` in the `cli_opts` module
    """
    if storage_policy_type != StoragePolicyType.CPG:
        raise PipelineError(f'Unsupported storage policy {storage_policy_type}')

    storage_provider = CpgStorageProvider(cloud)
    cohort = Cohort(
        analysis_dataset_name=analysis_dataset,
        namespace=namespace,
        name=name,
        storage_provider=storage_provider,
    )
    refs = RefData(storage_provider.get_ref_bucket())

    status_reporter: StatusReporter | None = None
    input_provider: InputProvider | None = None
    if (
        input_provider_type == InputProviderType.SMDB
        or status_reporter_type == StatusReporterType.SMDB
    ):
        smdb = SMDB(cohort.analysis_dataset.name)
        if status_reporter_type == StatusReporterType.SMDB:
            status_reporter = SmdbStatusReporter(smdb)
        if input_provider_type == InputProviderType.SMDB:
            input_provider = SmdbInputProvider(smdb)

    if input_provider_type == InputProviderType.CSV:
        if not input_csv:
            raise PipelineError(
                f'input_csv (--input-csv) should be provided '
                f'with input_provider_type=InputProviderType.CSV '
                f'(--input-provider {InputProviderType.CSV.value})'
            )
        input_provider = CsvInputProvider(to_path(input_csv).open())

    if not description:
        description = name
        if version:
            description += f' {version}'
        if datasets:
            datasets_ = set(datasets)
            if skip_datasets:
                datasets_ -= set(skip_datasets or [])
            description += ': ' + ', '.join(datasets_)

    tmp_bucket = to_path(
        cohort.analysis_dataset.get_tmp_bucket(
            version=(name + (f'/{version}' if version else ''))
        )
    )
    hail_billing_project = get_billing_project(cohort.analysis_dataset.stack)
    hail_bucket = tmp_bucket / 'hail'
    batch: Batch = setup_batch(
        description=description,
        billing_project=hail_billing_project,
        hail_bucket=hail_bucket,
    )

    if datasets and input_provider:
        input_provider.populate_cohort(
            cohort=cohort,
            dataset_names=datasets,
            skip_samples=skip_samples,
            only_samples=only_samples,
            skip_datasets=skip_datasets,
            ped_files=ped_files,
        )

    if force_samples:
        for s in cohort.get_samples():
            if s.id in force_samples:
                logger.info(f'Force rerunning sample {s.id} even if its outputs exist')
                s.forced = True

    return Pipeline(
        cohort=cohort,
        namespace=namespace,
        reference_data=refs,
        batch=batch,
        hail_billing_project=hail_billing_project,
        hail_bucket=hail_bucket,
        status_reporter=status_reporter,
        stages=stages,
        first_stage=first_stage,
        last_stage=last_stage,
        version=version,
        check_intermediates=check_intermediates,
        check_expected_outputs=check_expected_outputs,
        skip_samples_with_missing_input=skip_samples_with_missing_input,
        config=config,
        tmp_bucket=tmp_bucket,
        local_dir=local_dir,
        dry_run=dry_run,
        keep_scratch=keep_scratch,
    )
