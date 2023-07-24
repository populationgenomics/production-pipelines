from dataclasses import asdict, dataclass, field
from typing import Any, Literal, Optional

from cpg_utils import Path

from .dataset import DEFAULT_DATASET_NAME
from .types import DatasetId, SequencingGroupId, SequencingType, StageName


def _remove_none_values(d: dict[str, Any]) -> dict[str, Any]:
    """
    Recursively remove keys with None values from a dictionary.
    """
    return {
        key: (_remove_none_values(value) if isinstance(value, dict) else value)
        for key, value in d.items()
        if value is not None
    }


@dataclass(kw_only=True)
class WorkflowConfig:
    """
    See `<root>/cpg_workflows/defaults.toml` for a description of these options.

    Note: that the `Optional` type is used to indicate that a key is optional, but
    some stages may require it.
    """

    # ---- Core
    # The ID of this analysis', which may be an analysis over all datasets in the
    # `input_datasets` list spanning over different projects
    dataset: DatasetId = DEFAULT_DATASET_NAME

    # Usually 'genome' or 'exome', but will probably be expanded over time
    sequencing_type: SequencingType = 'genome'

    # GCP bucket access level. Use 'test' for all tests
    access_level: Literal['main', 'test'] = 'test'

    # Possibly unused, or used by cpg_utils somewhere?
    dataset_gcp_project: Optional[str] = None

    # Description of the workflow (to display in the Batch GUI)
    description: Optional[str] = None

    # Name of the workflow passed to Hail Batch instance
    name: Optional[str] = None

    # List of dataset IDs that will be processed during this analysis. Leave blank
    # if the input dataset is the dataset being processed.
    input_datasets: Optional[list[DatasetId]] = None

    # Version of this analysis... Not sure versioning is implemented yet.
    output_version: Optional[str] = None

    # Reference fasta file path to use instead of Broad's reference.
    ref_fasta: Optional[str | Path] = None

    # Skip quality control stages
    skip_qc: Optional[bool] = None

    # We only support one status reporter right now, and that is Metamist.
    status_reporter: Optional[Literal['metamist']] = None

    path_scheme: Optional[str] = None  # Possibly unused
    local_dir: Optional[str] = None  # Possibly unused

    # ---- Stage specific resources
    resources: Optional[dict[StageName, dict[str, Any]]] = None

    # ---- Technical
    # The driver image use when running an analysis on Hail Batch. Not used locally.
    driver_image: Optional[str] = None

    # Only computes DAG of stages to run, but doesn't call methods to queue jobs.
    dry_run: Optional[bool] = None

    # Request high memory compute instances
    highmem_workers: Optional[bool] = None

    # ---- Checks
    # Check input file existence (e.g. FASTQ files). When they are missing, the
    # `skip_sgs_with_missing_input` option controls whether such sequencing groups
    #  should be ignored, or it should cause raising an error. Defaults to `True`.
    check_inputs: Optional[bool] = None

    # Within jobs, check all in-job intermediate files for possible reuse. If set to
    # `False`, will overwrite all intermediates. Used by `utils.can_reuse(path)`.
    # Defaults to `True` in `utils.can_reuse(path)`, but stages may override this.
    check_intermediates: Optional[bool] = None

    # Before running a stage, check if its outputs already exist. If they exist,
    # then do not run this stage. Defaults to `False`.
    check_expected_outputs: Optional[bool] = None

    # ---- Stages
    first_stages: Optional[list[StageName]] = None
    last_stages: Optional[list[StageName]] = None
    only_stages: Optional[list[StageName]] = None
    skip_stages: Optional[list[StageName]] = None
    force_stages: Optional[list[StageName]] = None
    # pylint: disable=invalid-name
    allow_missing_outputs_for_stages: Optional[list[StageName]] = None

    # ---- Dataset options
    skip_datasets: Optional[list[DatasetId]] = None

    # ---- Sequencing group options
    skip_sgs: Optional[list[SequencingGroupId]] = None
    only_sgs: Optional[list[SequencingGroupId]] = None
    force_sgs: Optional[list[SequencingGroupId]] = None
    skip_stages_for_sgs: Optional[dict[StageName, list[SequencingGroupId]]] = None

    # For the first (not-skipped) stage, if the input for a sequencing group does not
    # exist, skip this instance instead of failing. For example, if the first stage is
    # Align, and `sequencing_group.alignment_input` for a sequencing group does not
    # exist, remove this sequencing group (i.e the `active` attiribute will be set to
    # `False`), instead of failing. In other words, ignore sequencing groups that are
    # missing results from skipped stages that the non-skipped stage might require.
    skip_sgs_with_missing_input: Optional[bool] = None

    # ---- Other
    realign_from_cram_version: Optional[str] = None
    cram_version_reference: Optional[dict[str, str | Path]] = None
    intervals_path: Optional[str | Path] = None
    reblock_gq_bands: Optional[list[int]] = None
    create_es_index_for_datasets: Optional[list[str]] = None
    scatter_count: Optional[int] = None
    scatter_count_genotype: Optional[int] = None
    vds_version: Optional[str] = None
    use_gnarly: Optional[bool] = None
    use_as_vqsr: Optional[bool] = None
    write_vcf: Optional[list[str]] = None


@dataclass(kw_only=True)
class HailConfig:
    """
    Configuration options used to configure Hail batch jobs.

    Attributes:
        dry_run (bool, optional):
            Runs Hail Batch in dry_run mode, meaning jobs are not executed. Defaults
            to `True`.

        backend (Optional[Literal['batch', 'local']], optional):
            Specifies which backend Hail Batch is initialised with. Use 'local' for
            local testing, and 'batch' for running on Google Cloud. Defaults to `local`.

        query_backend (Literal['spark', 'batch', 'local', 'spark_local'], optional):
            Specifies which backend Hail is initialised with when calling
            `cpg_utils.hail_batch.start_query_context`. Defaults to `spark`.

        billing_project (Optional[str], optional):
            The GCP billing project that Hail Batch will use to provision resources,
            and run jobs. Required if `backend` is set to `batch`. Defaults to `None`.

        pool_label (Optional[str], optional):
            Sets preemptible Hail Batch jobs' `_pool_label` attribute to this value.
            Defaults to `None`.

        delete_scratch_on_exit (Optional[bool], optional):
            If True, delete temporary directories containing intermediate files after
            the batch has finished executing. Defaults to `None`.

        cancel_after_n_failures (Optional[int], optional):
            Automatically cancel the batch after N failures have occurred. The default
            behavior is there is no limit on the number of failures. Only applicable
            if `backend` is set to `batch`. Must be greater than 0. Defaults to `None`.

        default_memory (Optional[str], optional):
            Memory setting to use by default if not specified by a Hail Batch Job
            instance. Defaults to `None`.

        default_timeout (Optional[str], optional):
            Maximum time in seconds for a job to run before being killed. Only
            applicable if `backend` is set to `batch`. There is no timeout if this
            value is not set. Defaults to `None`.
    """

    dry_run: bool = True
    backend: Literal['batch', 'local'] = 'local'
    query_backend: Literal['spark', 'batch', 'local', 'spark_local'] = 'spark'
    billing_project: Optional[str] = None
    pool_label: Optional[str] = None
    delete_scratch_on_exit: Optional[bool] = True
    cancel_after_n_failures: Optional[int] = None
    default_memory: Optional[str] = None
    default_timeout: Optional[str] = None


@dataclass(kw_only=True)
class StorageConfig:
    """
    Technically `cpg_utils` can pull from any storge key in a workflow config, but
    these are the keys that we generally stick with. For example:

    [storage.fewgenomes]
    default = "gs://cpg-fewgenomes-main"
    analysis = "gs://cpg-fewgenomes-main-analysis"
    tmp = "gs://cpg-fewgenomes-main-tmp"
    web = "gs://cpg-fewgenomes-main-web"
    web_url = "https://main-web.populationgenomics.org.au/fewgenomes"
    """

    default: Optional[str | Path] = None
    analysis: Optional[str | Path] = None
    tmp: Optional[str | Path] = None
    web: Optional[str | Path] = None
    web_url: Optional[str | Path] = None


def default_hail_config() -> HailConfig:
    # Wrapped in function to avoid accidental mutation of default values
    return HailConfig(
        dry_run=True,
        backend='local',
        query_backend='spark',
        delete_scratch_on_exit=True,
    )


def default_workflow_config() -> WorkflowConfig:
    # Wrapped in function to avoid accidental mutation of default values
    return WorkflowConfig(
        dataset=DEFAULT_DATASET_NAME,
        access_level='test',
        sequencing_type='genome',
        check_inputs=True,
    )


@dataclass(kw_only=True)
class PipelineConfig:
    """
    Utility to create either a config `dict` or TOML string. There is subect to change
    as configration is refactored. There is an `**other` attribute to allow a catch-all
    for any configuration section which are not present heres explicit keyword
    arguments. The config keys not represented here that should go in `other` are:

        - references
            - broad
            - gnomad
            - gatk_sv
        - validation.sample_map
        - combiner
        - large_cohort
        - sv_ref_panel
        - vqsr
        - cramqc
        - qc_thresholds
        - elasticsearch
        - slack
        - stripy

    Lastly, if `keep_dict_keys_with_none` and `as_dict` are `True`, keys which map
    to `None` will be kept in the output `dict`, otherwise they are removed. Note that
    TOML removes keys with `None` values by default, so `keep_dict_keys_with_none`
    will have no effect if `as_dict` is `False`.

    Attributes:
        workflow (WorkflowConfig, optional):
            A `WorkflowConfig` dictionary instance. Defaults to the a minimal
            configuration for local testing.

        hail (HailConfig, optional):
            A `HailConfig` dictionary instance. Defaults to a config to enable local
            testing.

        images: (dict[str, str | Path], optional):
            A dictionary of docker image names and their corresponding image paths.
            Defaults to `{}`.

        storage (dict[DatasetId, StorageConfig], optional):
            Storage paths for a dataset. Each dataset's storage configuration can
            theoretically accept any key, but we typically use 'default', 'tmp', 'web',
            'analysis', 'web_url' paths. Defaults to `{}`.

        other (dict[str, dict[str, Any]], optional):
            A catch-all for any configuration section which are not present here as
            explicit keyword arguments. Defaults to `{}`.
    """

    workflow: WorkflowConfig = field(default_factory=default_workflow_config)
    hail: HailConfig = field(default_factory=default_hail_config)
    images: dict[str, str | Path] = field(default_factory=dict)
    storage: dict[DatasetId, StorageConfig] = field(default_factory=dict)
    other: dict[str, dict[str, Any]] = field(default_factory=dict)

    def set_storage(
        self, dataset: DatasetId, storage: StorageConfig
    ) -> 'PipelineConfig':
        self.storage[dataset] = storage
        return self

    def set_image(self, name: str, path: str | Path) -> 'PipelineConfig':
        self.images[name] = path
        return self

    def set_other(self, key: str, value: dict[str, Any]) -> 'PipelineConfig':
        self.other[key] = value
        return self

    def as_dict(self) -> dict[str, Any]:
        d = _remove_none_values(asdict(self))

        # Splat `other` so that these keys are accessible at the root level of the dict
        # required for the TOML config to work properly.
        other = d.pop('other')

        # Leave these out unless they contain something.
        storage = d.pop('storage')
        images = d.pop('images')

        config = {**d, **other}

        if storage:
            config['storage'] = storage
        if images:
            config['images'] = images

        return config
