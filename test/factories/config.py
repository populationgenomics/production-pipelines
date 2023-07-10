from dataclasses import asdict, dataclass
from typing import Any, Optional

import toml
from cpg_utils import Path

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


class TomlAnyPathEncoder(toml.TomlEncoder):
    """Support for CPG path objects in TOML"""

    def dump_value(self, v):
        if isinstance(v, Path):
            v = str(v)
        return super().dump_value(v)


@dataclass
class DictionaryMixin:
    """Mixin to convert dataclass to dictionary."""

    def as_dict(self, keep_keys_with_none: bool = False) -> dict:
        """
        Return dataclass as a dictionary.

        Args:
            keep_keys_with_none (bool, optional): Remove keys which map to `None`.
            Defaults to False.

        Returns:
            dict[str, Any]
        """
        config = asdict(self)
        if not keep_keys_with_none:
            config = _remove_none_values(config)
        return config


@dataclass(kw_only=True)
class WorkflowConfig(DictionaryMixin):
    """
    See `<root>/cpg_workflows/defaults.toml` for a description of these options.
    """

    # ---- Core
    dataset: DatasetId = 'test'
    access_level: str = 'test'
    sequencing_type: SequencingType = 'genome'
    # Possibly unused, or used by cpg_utils somewhere?
    dataset_gcp_project: Optional[str] = None
    # Description of the workflow (to display in the Batch GUI)
    description: Optional[str] = None
    # Name of the workflow passed to Hail Batch instance
    name: Optional[str] = None
    input_datasets: Optional[list[DatasetId]] = None
    output_version: Optional[str] = None
    ref_fasta: Optional[str | Path] = None
    skip_qc: Optional[bool] = None
    status_reporter: Optional[str] = None
    path_scheme: Optional[str] = None  # Possibly unused
    local_dir: Optional[str] = None  # Possibly unused

    # ---- Stage specific resources
    resources: Optional[dict[StageName, dict[str, Any]]] = None

    # ---- Technical
    driver_image: Optional[str] = None
    dry_run: Optional[bool] = None
    highmem_workers: Optional[bool] = None

    # ---- Checks
    # Check input file existence (e.g. FASTQ files). When they are missing,
    # the `skip_sgs_with_missing_input` option controls whether such
    # sequencing groups should be ignored, or it should cause raising an error.
    check_inputs: Optional[bool] = None
    # Within jobs, check all in-job intermediate files for possible reuse.
    # If set to False, will overwrite all intermediates. Used by `utils.can_reuse(path)`.
    check_intermediates: Optional[bool] = None
    # Before running a stage, check if its inputs (i.e.the expected outputs from
    # all required stages) already exist. If it exists, do not submit stage jobs.
    check_expected_outputs: Optional[bool] = None

    # ---- Stages
    first_stages: Optional[list[StageName]] = None
    last_stages: Optional[list[StageName]] = None
    only_stages: Optional[list[StageName]] = None
    skip_stages: Optional[list[StageName]] = None
    force_stages: Optional[list[StageName]] = None
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
class HailConfig(DictionaryMixin):
    dry_run: bool = True
    backend: str = 'local'
    query_backend: Optional[str] = None
    billing_project: Optional[str] = None
    pool_label: Optional[str] = None
    delete_scratch_on_exit: Optional[bool] = None


@dataclass(kw_only=True)
class StorageConfig(DictionaryMixin):
    default: Optional[str | Path] = None
    analysis: Optional[str | Path] = None
    tmp: Optional[str | Path] = None
    web: Optional[str | Path] = None
    web_url: Optional[str | Path] = None


def create_config(
    workflow: WorkflowConfig = WorkflowConfig(),  # noqa: B008
    hail: HailConfig = HailConfig(),  # noqa: B008
    images: Optional[dict[str, str]] = None,
    storage: Optional[dict[str, StorageConfig]] = None,
    as_dict: bool = False,
    keep_dict_keys_with_none: bool = False,
    **other,
) -> str | dict[str, Any]:
    """
    Utility to create either a config `dict` or TOML string. There is subect to change
    as configration is refactored. There is a `**other` argument to allow for a
    catch-all for configuration fields which are not present here s explicit keyword
    arguments. If `keep_dict_keys_with_none` and `as_dict` are `True`, keys which map
    to `None` will be kept in the output dictionary.

    Returns:
        str | dict[str, Any]
    """

    # Inform user of incorrect usage of `keep_dict_keys_with_none`, since TOML will
    # remove keys which map to `None` by default, and I don't want to override that
    # behaviour right now.
    if keep_dict_keys_with_none and not as_dict:
        raise ValueError(
            "'keep_dict_keys_with_none' can only be set to 'True' if "
            + "'as_dict' is also 'True'"
        )

    config = {
        'workflow': workflow.as_dict(keep_dict_keys_with_none),
        'hail': hail.as_dict(keep_dict_keys_with_none),
        'images': images,
        'storage': storage,
        **other,
    }

    # Treating None as distinct from an empty dict because we want to be able to
    # test what happens if keys are missing from a potentially empty section.
    if not keep_dict_keys_with_none:
        config = _remove_none_values(config)

    if as_dict:
        return config

    return toml.dumps(config, encoder=TomlAnyPathEncoder())
