from dataclasses import asdict, dataclass, field
from typing import Any, Literal

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

    Note: that the `Union` type with `None` is used to indicate that a key is
    optional, but some stages may require it.

    If you find a config key not present here, please add it along with some
    documnentation on what it is, and what it is used for.
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
    dataset_gcp_project: str | None = None

    # Description of the workflow (to display in the Batch GUI)
    description: str | None = None

    # Name of the workflow passed to Hail Batch instance
    name: str | None = None

    # List of dataset IDs that will be processed during this analysis. Leave blank
    # if the input dataset is the dataset being processed.
    input_datasets: list[DatasetId] | None = None

    # List of cohort IDs that will be processed during this analysis.
    input_cohorts: list[str] | None = None

    # Version of this analysis... Not sure versioning is implemented yet.
    output_version: str | None = None

    # Reference fasta file path to use instead of Broad's reference.
    ref_fasta: str | Path | None = None

    # Skip quality control stages
    skip_qc: bool | None = None

    # We only support one status reporter right now, and that is Metamist.
    status_reporter: Literal['metamist'] | None = None

    path_scheme: str | None = None  # Possibly unused
    local_dir: str | None = None  # Possibly unused

    # ---- Stage specific resources
    resources: dict[StageName, dict[str, Any]] | None = None

    # ---- Technical
    # The driver image use when running an analysis on Hail Batch. Not used locally.
    driver_image: str | None = None

    # Only computes DAG of stages to run, but doesn't call methods to queue jobs.
    dry_run: bool | None = None

    # Request high memory compute instances
    highmem_workers: bool | None = None

    # ---- Checks
    # Check input file existence (e.g. FASTQ files). When they are missing, the
    # `skip_sgs_with_missing_input` option controls whether such sequencing groups
    #  should be ignored, or it should cause raising an error. Defaults to `True`.
    check_inputs: bool | None = None

    # Within jobs, check all in-job intermediate files for possible reuse. If set to
    # `False`, will overwrite all intermediates. Used by `utils.can_reuse(path)`.
    # Defaults to `True` in `utils.can_reuse(path)`, but stages may override this.
    check_intermediates: bool | None = None

    # Before running a stage, check if its outputs already exist. If they exist,
    # then do not run this stage. Defaults to `False`.
    check_expected_outputs: bool | None = None

    # ---- Stages
    first_stages: list[StageName] | None = None
    last_stages: list[StageName] | None = None
    only_stages: list[StageName] | None = None
    skip_stages: list[StageName] | None = None
    force_stages: list[StageName] | None = None
    # pylint: disable=invalid-name
    allow_missing_outputs_for_stages: list[StageName] | None = None

    # ---- Dataset options
    skip_datasets: list[DatasetId] | None = None

    # ---- Sequencing group options
    skip_sgs: list[SequencingGroupId] | None = None
    only_sgs: list[SequencingGroupId] | None = None
    force_sgs: list[SequencingGroupId] | None = None
    skip_stages_for_sgs: dict[StageName, list[SequencingGroupId]] | None = None

    # For the first (not-skipped) stage, if the input for a sequencing group does not
    # exist, skip this instance instead of failing. For example, if the first stage is
    # Align, and `sequencing_group.alignment_input` for a sequencing group does not
    # exist, remove this sequencing group (i.e the `active` attiribute will be set to
    # `False`), instead of failing. In other words, ignore sequencing groups that are
    # missing results from skipped stages that the non-skipped stage might require.
    skip_sgs_with_missing_input: bool | None = None

    # ---- Other
    realign_from_cram_version: str | None = None
    cram_version_reference: dict[str, str | Path] | None = None
    intervals_path: str | Path | None = None
    reblock_gq_bands: list[int] | None = None
    create_es_index_for_datasets: list[str] | None = None
    scatter_count: int | None = None
    scatter_count_genotype: int | None = None
    vds_version: str | None = None
    use_gnarly: bool | None = None
    use_as_vqsr: bool | None = None
    write_vcf: list[str] | None = None


@dataclass(kw_only=True)
class HailConfig:
    """
    Configuration options used to configure Hail batch jobs.

    Attributes:
        dry_run (bool, optional):
            Runs Hail Batch in dry_run mode, meaning jobs are not executed. Defaults
            to `True`.

        backend (Literal['batch', 'local' | None], optional):
            Specifies which backend Hail Batch is initialised with. Use 'local' for
            local testing, and 'batch' for running on Google Cloud. Defaults to `local`.

        query_backend (Literal['spark', 'batch', 'local', 'spark_local'], optional):
            Specifies which backend Hail is initialised with when calling
            `cpg_utils.hail_batch.start_query_context`. Defaults to `spark_local`.

        billing_project (str | None, optional):
            The GCP billing project that Hail Batch will use to provision resources,
            and run jobs. Required if `backend` is set to `batch`. Defaults to `None`.

        pool_label (str | None, optional):
            Sets preemptible Hail Batch jobs' `_pool_label` attribute to this value.
            Defaults to `None`.

        delete_scratch_on_exit (bool | None, optional):
            If True, delete temporary directories containing intermediate files after
            the batch has finished executing. Defaults to `None`.

        cancel_after_n_failures (int | None, optional):
            Automatically cancel the batch after N failures have occurred. The default
            behavior is there is no limit on the number of failures. Only applicable
            if `backend` is set to `batch`. Must be greater than 0. Defaults to `None`.

        default_memory (str | None, optional):
            Memory setting to use by default if not specified by a Hail Batch Job
            instance. Defaults to `None`.

        default_timeout (str | None, optional):
            Maximum time in seconds for a job to run before being killed. Only
            applicable if `backend` is set to `batch`. There is no timeout if this
            value is not set. Defaults to `None`.
    """

    dry_run: bool = True
    backend: Literal['batch', 'local'] = 'local'
    query_backend: Literal['spark', 'batch', 'local', 'spark_local'] = 'spark_local'
    billing_project: str | None = None
    pool_label: str | None = None
    delete_scratch_on_exit: bool | None = True
    cancel_after_n_failures: int | None = None
    default_memory: str | None = None
    default_timeout: str | None = None


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

    default: str | Path | None = None
    analysis: str | Path | None = None
    tmp: str | Path | None = None
    web: str | Path | None = None
    web_url: str | Path | None = None


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

        - validation.sample_map
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

        storage (dict[DatasetId, StorageConfig], optional):
            Storage paths for a dataset. Each dataset's storage configuration can
            theoretically accept any key, but we typically use 'default', 'tmp', 'web',
            'analysis', 'web_url' paths. Defaults to `{}`.

        images: (dict[str, str | Path], optional):
            A dictionary of docker image names and their corresponding image paths.
            Defaults to `{}`.

        references (dict[str, str | dict[str, Any]], optional):
            A dictionary containing reference options. May also contain nested
            dictionaries for broad, gnomad, gatk_sv and other references. Defaults to
            `{}`.

        large_cohort (dict[str, Any], optional):
            A dictionary containing large cohort pipeline configuration parameters and
            values used by functions in the large cohort module (e.g. `n_pcs`,
            `sample_qc_cutoffs`, `combiner`). Defaults to `{}`.

        other (dict[str, dict[str, Any]], optional):
            A catch-all for any configuration section which are not present here as
            explicit keyword arguments. Defaults to `{}`.
    """

    workflow: WorkflowConfig = field(default_factory=default_workflow_config)
    hail: HailConfig = field(default_factory=default_hail_config)
    storage: dict[DatasetId, StorageConfig] = field(default_factory=dict)
    images: dict[str, str | Path] = field(default_factory=dict)
    references: dict[str, Any | dict[str, Any]] = field(default_factory=dict)
    large_cohort: dict[str, Any] = field(default_factory=dict)
    other: dict[str, dict[str, Any]] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        return self.as_dict()[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.as_dict().get(key, default)

    def set_storage(self, dataset: DatasetId, storage: StorageConfig) -> 'PipelineConfig':
        self.storage[dataset] = storage
        return self

    def set_image(self, name: str, path: str | Path) -> 'PipelineConfig':
        self.images[name] = path
        return self

    def set_references(self, name: str, value: str | dict[str, Any]) -> 'PipelineConfig':
        self.references[name] = value
        return self

    def set_large_cohort(self, key: str, value: Any) -> 'PipelineConfig':
        self.large_cohort[key] = value
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
        references = d.pop('references')
        large_cohort = d.pop('large_cohort')

        config = {**d, **other}

        if storage:
            config['storage'] = storage
        if images:
            config['images'] = images
        if references:
            config['references'] = references
        if large_cohort:
            config['large_cohort'] = large_cohort

        return config
