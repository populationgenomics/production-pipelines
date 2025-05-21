"""
Runs the hail combiner
This is copied from the large_cohort implementation, as RD needs to alter some defaults
Specifically we need to drop/modify the Force argument to resume from previous Combiner plan/runs
"""


def run(
    output_vds_path: str,
    sequencing_type: str,
    tmp_prefix: str,
    genome_build: str,
    save_path: str | None,
    gvcf_paths: list[str] | None = None,
    vds_path: str | None = None,
    specific_intervals: list[str] | None = None,
    force_new_combiner: bool = False,
    sgs_to_remove: list[str] | None = None,
) -> None:
    """
    Runs the combiner -

    1. this method receives an existing VDS path or None - if None, combine from scratch. If populated, the combiner
       run will use that as a base, and add additional samples to it
    2. if this method receives sgs_to_remove and a vds_path, the first step is to create a new VDS from temp with the
       samples removed. This is then used as the base for the combiner run
    3. if this method receives gvcf_paths, the combiner will run with those as the input paths to add
    4. if there are no gVCFs, the VDS with samples removed is written to the output path directly

    Args:
        output_vds_path (str): eventual output path for the VDS
        sequencing_type (str): genome/exome, relevant in selecting defaults
        tmp_prefix (str): where to store temporary combiner intermediates
        genome_build (str): GRCh38
        save_path (str | None): where to store the combiner plan, or where to resume from
        gvcf_paths (list[str] | None): list of paths to GVCFs
        vds_path (str | None): a single VDS, or None - this is where a combiner can continue on from
        specific_intervals (list[str] | None): list of intervals to use for the combiner, if using non-standard
        force_new_combiner (bool): whether to force a new combiner run, or permit resume from a previous one
        sgs_to_remove (list[str] | None): list of sample groups to remove from the combiner
    """
    import logging

    import hail as hl

    from cpg_utils.config import config_retrieve
    from cpg_utils.hail_batch import init_batch
    from cpg_workflows.utils import exists

    # set up a quick logger inside the job
    logging.basicConfig(level=logging.INFO)

    init_batch(
        worker_memory=config_retrieve(['combiner', 'worker_memory']),
        driver_memory=config_retrieve(['combiner', 'driver_memory']),
        driver_cores=config_retrieve(['combiner', 'driver_cores']),
    )

    # Load from save, if supplied (log correctly depending on force_new_combiner)
    if save_path and force_new_combiner:
        logging.info(f'Combiner plan {save_path} will be ignored/written new')

    elif save_path:
        logging.info(f'Resuming combiner plan from {save_path}')

    if specific_intervals:
        logging.info(f'Using specific intervals: {specific_intervals}')
        intervals = hl.eval(
            [hl.parse_locus_interval(interval, reference_genome=genome_build) for interval in specific_intervals],
        )
    else:
        intervals = None

    if not (sgs_to_remove or gvcf_paths):
        raise ValueError('No samples to remove or gVCFs to add - please provide at least one of these')

    # logical steps -
    # 1. if there are samples to remove, do that first
    #   a. if there are no samples to add, just remove the samples and write to eventual output path
    #   b. if there are samples to add, write this to a temporary path, set force to true, and then run the combiner
    # 2. if there are samples to add, run the combiner with the final path as the output path

    # 1 - do we need to do removal?
    if vds_path and sgs_to_remove:
        logging.info(f'Removing sample groups {sgs_to_remove} from {vds_path}')

        temp_path = f'{tmp_prefix}/combiner_removal_temp.vds'

        # this will only exist if the previous removal was successful, AND we have additional gVCFs to add
        if exists(temp_path):
            logging.info(f'Found existing VDS at {temp_path}, skipping removal step')
            vds_path = temp_path

        else:
            vds = hl.vds.read_vds(vds_path)
            vds = hl.vds.filter_samples(vds, samples=sgs_to_remove, keep=False, remove_dead_alleles=True)

            # 1a. if there are no samples to add, just remove the samples and write to eventual output path
            if not gvcf_paths:
                logging.info(f'Writing to {output_vds_path}')
                vds.write(output_vds_path)
                return

            # 1b. if there are samples to add, write this to a temporary path, set force to true, then run the combiner
            else:
                # we've changed the VDS base, so we need to force a new combiner plan
                force_new_combiner = True

                logging.info(f'Writing with removed SGs to {temp_path}')
                vds.write(temp_path)
                vds_path = temp_path

    # 2 - do we need to run the combiner?
    combiner = hl.vds.new_combiner(
        output_path=output_vds_path,
        save_path=save_path,
        gvcf_paths=gvcf_paths,
        vds_paths=[vds_path] if vds_path else None,
        reference_genome=genome_build,
        temp_path=tmp_prefix,
        use_exome_default_intervals=sequencing_type == 'exome',
        use_genome_default_intervals=sequencing_type == 'genome',
        intervals=intervals,
        force=force_new_combiner,
        branch_factor=config_retrieve(['combiner', 'branch_factor']),
        target_records=config_retrieve(['combiner', 'target_records']),
        gvcf_batch_size=config_retrieve(['combiner', 'gvcf_batch_size']),
    )

    combiner.run()
