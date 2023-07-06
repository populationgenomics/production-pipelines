from typing import Any

import toml

from .types import SequencingType


def create_config(
    dataset: str = "test",
    dataset_gcp_project: str = "local-test",
    access_level: str = "test",
    sequencing_type: SequencingType = "genome",
    check_inputs: bool = False,
    cram_version_reference: dict[str, str] | None = None,
    realign_from_cram_version: str = "",
    workflow_resources: dict[str, Any] | None = None,
    images: dict[str, str] | None = None,
    references: dict[str, str | dict[str, str]] | None = None,
    resource_overrides: dict[str, Any] | None = None,
    storage: dict[str, dict[str, str]] | None = None,
    **other,
) -> str:
    config = f"""
    [workflow]
    dataset_gcp_project = "{dataset_gcp_project}"
    access_level = "{access_level}"
    dataset = "{dataset}"
    sequencing_type = "{sequencing_type}"
    realign_from_cram_version = "{realign_from_cram_version}"
    check_inputs = {str(check_inputs).lower()}

    [hail]
    dry_run = true
    backend = "local"

    [images]
    :images

    :references

    [cram_version_reference]
    :cram_ref_map

    [resource_overrides]
    :resource_overrides
    
    [workflow.resources.Align]
    :align_resources

    :storage

    :other
    """
    if not references:
        references = {
            "broad": {
                "dragmap_prefix": "gs://a-cpg-bucket/dragen_reference/",
                "ref_fasta": "gs://a-cpg-bucket/hg38.fasta",
            }
        }

    if not images:
        images = {
            "dragmap": "dragmap_image:1.3.0",
            "biobambam2": "biobambam2_image:2.0.87",
            "picard": "picard_image:2.27.4",
            "samtools": "samtools_image:1.11",
        }

    config = config.replace(":cram_ref_map", toml.dumps(cram_version_reference or {}))
    config = config.replace(":align_resources", toml.dumps(workflow_resources or {}))
    config = config.replace(":references", toml.dumps({"references": references}))
    config = config.replace(":images", toml.dumps(images))
    config = config.replace(":storage", toml.dumps(storage or {}))
    config = config.replace(":resource_overrides", toml.dumps(resource_overrides or {}))

    # Catch all for all other config sections
    config = config.replace(":other", toml.dumps(other))

    return config
