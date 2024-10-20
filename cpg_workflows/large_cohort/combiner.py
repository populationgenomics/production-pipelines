"""Runs the combiner."""

from typing import TYPE_CHECKING

from hail.vds import new_combiner

from cpg_utils import Path
from cpg_utils.hail_batch import init_batch

if TYPE_CHECKING:  # TCH002 https://docs.astral.sh/ruff/rules/typing-only-third-party-import/
    from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner


def run_combiner(
    output_vds_path: Path,
    sequencing_type: str,
    tmp_prfx: str,
    gvcf_paths: list[str] | None = None,
    vds_paths: list[str] | None = None,
) -> None:
    use_genome_default_intervals: bool = False
    use_exome_default_intervals: bool = False

    init_batch()

    if sequencing_type == "genome":
        use_genome_default_intervals = True
    elif sequencing_type == "exome":
        use_exome_default_intervals = True

    combiner: VariantDatasetCombiner = new_combiner(
        output_path=str(output_vds_path),
        gvcf_paths=gvcf_paths,
        vds_paths=vds_paths,
        reference_genome="GRCh38",
        temp_path=tmp_prfx,
        use_exome_default_intervals=use_exome_default_intervals,
        use_genome_default_intervals=use_genome_default_intervals,
        force=True,
    )

    combiner.run()
