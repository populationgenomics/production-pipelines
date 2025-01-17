"""Runs the hail combiner

Inputs:
    - output_vds_path: str - The destination that the VDS will be saved at
    - sequencing_type: str - Used to specify what intervals to use (exome or genome)
    - tmp_prefix: str - Where to store temporary combiner plans
    - genome_build: str - What reference genome to use for the combiner
    - gvcf_paths: list[str] | None - The optional list of gvcf paths in string format to combine
    - vds_paths: list[str] | None - The optional list of VDS paths in string format to combine
"""

from typing import TYPE_CHECKING

from hail.vds import new_combiner
from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import can_reuse, to_path


def run(
    output_vds_path: str,
    sequencing_type: str,
    tmp_prefix: str,
    genome_build: str,
    save_path: str | None,
    force_new_combiner: bool = False,
    gvcf_paths: list[str] | None = None,
    vds_paths: list[str] | None = None,
    specific_intervals: list[str] | None = None,
) -> None:
    """
    Runs the combiner

    Args:
        output_vds_path (str): eventual output path for the VDS
        sequencing_type (str): genome/exome, relevant in selecting defaults
        tmp_prefix (str): where to store temporary combiner intermediates
        genome_build (str): GRCh38
        save_path (str | None): where to store the combiner plan, or where to resume from
        force_new_combiner (bool): whether to force a new combiner run, or permit resume from a previous one
        gvcf_paths (list[str] | None): list of paths to GVCFs
        vds_paths (list[str] | None): list of paths to VDSs
        specific_intervals (list[str] | None): list of intervals to use for the combiner, if using non-standard
    """
    import logging

    import hail as hl

    # set up a quick logger inside the job
    logging.basicConfig(level=logging.INFO)

    if not can_reuse(to_path(output_vds_path)):
        init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)

        if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
            override_jar_spec(jar_spec)

        # Load from save, if supplied
        if save_path:
            if force_new_combiner:
                logging.info(f'Combiner plan {save_path} will be ignored/written new')
            else:
                logging.info(f'Resuming combiner plan from {save_path}')

        if specific_intervals:
            logging.info(f'Using specific intervals: {specific_intervals}')

            intervals = hl.eval(
                [hl.parse_locus_interval(interval, reference_genome=genome_build) for interval in specific_intervals],
            )

        else:
            intervals = None

        combiner: VariantDatasetCombiner = new_combiner(
            output_path=output_vds_path,
            save_path=save_path,
            gvcf_paths=gvcf_paths,
            vds_paths=vds_paths,
            reference_genome=genome_build,
            temp_path=tmp_prefix,
            use_exome_default_intervals=sequencing_type == 'exome',
            use_genome_default_intervals=sequencing_type == 'genome',
            intervals=intervals,
            force=force_new_combiner,
            branch_factor=config_retrieve(
                ['combiner', 'branch_factor'],
                VariantDatasetCombiner._default_branch_factor,
            ),
            gvcf_batch_size=config_retrieve(
                ['combiner', 'gvcf_batch_size'],
                None,
            ),
        )

        combiner.run()
