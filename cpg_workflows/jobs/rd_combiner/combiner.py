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
    vds_paths: list[str] | None = None,
    specific_intervals: list[str] | None = None,
    force_new_combiner: bool = False,
) -> None:
    """
    Runs the combiner

    Args:
        output_vds_path (str): eventual output path for the VDS
        sequencing_type (str): genome/exome, relevant in selecting defaults
        tmp_prefix (str): where to store temporary combiner intermediates
        genome_build (str): GRCh38
        save_path (str | None): where to store the combiner plan, or where to resume from
        gvcf_paths (list[str] | None): list of paths to GVCFs
        vds_paths (list[str] | None): list of paths to VDSs
        specific_intervals (list[str] | None): list of intervals to use for the combiner, if using non-standard
        force_new_combiner (bool): whether to force a new combiner run, or permit resume from a previous one
    """
    import logging

    import hail as hl
    from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner

    from cpg_utils.config import config_retrieve
    from cpg_utils.hail_batch import init_batch
    from cpg_workflows.batch import override_jar_spec

    # set up a quick logger inside the job
    logging.basicConfig(level=logging.INFO)

    init_batch(
        worker_memory=config_retrieve(['combiner', 'worker_memory'], 'highmem'),
        driver_memory=config_retrieve(['combiner', 'driver_memory'], 'highmem'),
        driver_cores=config_retrieve(['combiner', 'driver_cores'], 2),
    )
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revisions', 'combiner'], False):
        override_jar_spec(jar_spec)

    # Load from save, if supplied (log correctly depending on force_new_combiner)
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

    combiner = hl.vds.new_combiner(
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
        # we're defaulting to the protected class attributes here, which looks like a hack...
        # for branch factor and target records, the argument uses a specific value as a default
        # so if we don't find an entry in config, we can't pass None to the constructor...
        # we either access the protected class attributes, hard-code the default on our side,
        # or have two separate constructors depending on whether we override the default or not
        branch_factor=config_retrieve(
            ['combiner', 'branch_factor'],
            VariantDatasetCombiner._default_branch_factor,
        ),
        target_records=config_retrieve(
            ['combiner', 'target_records'],
            VariantDatasetCombiner._default_target_records,
        ),
        # this argument does default to None, and will be set to the default values within the constructor
        # so we're happy to pass None, no need to access the protected class attributes
        gvcf_batch_size=config_retrieve(
            ['combiner', 'gvcf_batch_size'],
            None,
        ),
    )

    combiner.run()
