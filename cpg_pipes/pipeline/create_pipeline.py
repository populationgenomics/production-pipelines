"""
Provides a `create_pipeline` function that create an instance of a `Pipeline`
using concrete implementations of the input and storage providers.
"""
import logging
from typing import Callable

from .exceptions import PipelineError
from .pipeline import Pipeline, Stage
from .. import Path, to_path
from ..providers import (
    StoragePolicyType,
    StatusReporterType,
    InputProviderType,
)
from ..providers.storage import Namespace, Cloud
from ..providers.cpg import (
    CpgStorageProvider,
    CpgStatusReporter,
    SmdbInputProvider,
)
from ..providers.cpg.smdb import SMDB
from ..providers.inputs import InputProvider, CsvInputProvider
from ..providers.status import StatusReporter
from ..targets import Dataset

logger = logging.getLogger(__file__)


def init_storage_provider(
    storage_policy_type: StoragePolicyType = None,
    cloud: Cloud = Cloud.GS,
):
    """
    Create storage provider based on enum options
    """
    if storage_policy_type == StoragePolicyType.CPG:
        return CpgStorageProvider(cloud)
    else:
        raise PipelineError(f'Unsupported storage policy {storage_policy_type}')


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
    check_inputs: bool = False,
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
    local_dir: Path | None = None,
) -> 'Pipeline':
    """
    Create a Pipeline instance. All options correspond to command line parameters
    described in `pipeline_click_options` in the `cli_opts` module
    """
    storage_provider = init_storage_provider(storage_policy_type, cloud)

    status_reporter: StatusReporter | None = None
    input_provider: InputProvider | None = None
    if (
        input_provider_type == InputProviderType.SMDB
        or status_reporter_type == StatusReporterType.CPG
    ):
        sm_proj = Dataset(analysis_dataset, namespace=namespace).stack
        smdb = SMDB(sm_proj)
        if status_reporter_type == StatusReporterType.CPG:
            status_reporter = CpgStatusReporter(smdb)
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

    return Pipeline(
        namespace=namespace,
        name=name.replace(' ', '_'),
        description=description,
        analysis_dataset_name=analysis_dataset,
        storage_provider=storage_provider,
        status_reporter=status_reporter,
        input_provider=input_provider,
        datasets=datasets,
        skip_datasets=skip_datasets,
        skip_samples=skip_samples,
        only_samples=only_samples,
        force_samples=force_samples,
        stages=stages,
        first_stage=first_stage,
        last_stage=last_stage,
        version=version,
        check_inputs=check_inputs,
        check_intermediates=check_intermediates,
        check_expected_outputs=check_expected_outputs,
        skip_samples_with_missing_input=skip_samples_with_missing_input,
        config=config,
        local_dir=local_dir,
        dry_run=dry_run,
        keep_scratch=keep_scratch,
    )
