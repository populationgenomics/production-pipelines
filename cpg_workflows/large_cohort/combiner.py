"""Runs the hail combiner

Inputs:
    - output_vds_path: str - The destination that the VDS will be saved at
    - sequencing_type: str - Used to specify what intervals to use (exome or genome)
    - tmp_prefix: str - Where to store temporary combiner plans
    - genome_build: str - What reference genome to use for the combiner
    - gvcf_paths: list[str] | None - The optional list of gvcf paths in string format to combine
    - vds_paths: list[str] | None - The optional list of VDS paths in string format to combine
"""

from typing import TYPE_CHECKING, Any

from hail.vds import new_combiner
from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner
from hailtop.batch.job import PythonJob

from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import get_batch, init_batch
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.utils import can_reuse, slugify, to_path
from metamist.graphql import gql, query

if TYPE_CHECKING:
    from graphql import DocumentNode

    from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner
    from hail.vds.variant_dataset import VariantDataset

    import cpg_utils


def _initalise_combiner_job(cohort: Cohort) -> PythonJob:
    j: PythonJob = get_batch().new_python_job('Combiner', (cohort.get_job_attrs() or {}) | {'tool': 'hail query'})  # type: ignore[reportUnknownArgumentType]
    j.image(config_retrieve(['workflow', 'driver_image']))
    j.memory(config_retrieve(['combiner', 'memory']))
    j.storage(config_retrieve(['combiner', 'storage']))

    # set this job to be non-spot (i.e. non-preemptible)
    # previous issues with preemptible VMs led to multiple simultaneous QOB groups processing the same data
    j.spot(config_retrieve(['combiner', 'preemptible_vms'], False))
    return j


def get_vds_ids_output(vds_id: int) -> tuple[str, list[str]]:
    get_vds_analysis_query: DocumentNode = gql(
        """
            query getVDSByAnalysisIds($vds_id: Int!) {
                analyses(id: {eq: $vds_id}) {
                    output
                    sequencingGroups {
                        id
                    }
                }
            }
        """,
    )
    query_results: dict[str, Any] = query(get_vds_analysis_query, variables={'vds_id': vds_id})
    vds_path: str = query_results['analyses'][0]['output']
    vds_sgids: list[str] = [sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']]
    return (vds_path, vds_sgids)


def combiner(cohort: Cohort, output_vds_path: str, save_path: str) -> PythonJob:
    workflow_config: dict[str, Any] = config_retrieve('workflow')
    combiner_config: dict[str, Any] = config_retrieve('combiner')

    combiner_tmp_path: cpg_utils.Path = (
        cohort.analysis_dataset.tmp_prefix()
        / 'vds'
        / f'{cohort.name}'
        / f'{cohort.id}-{combiner_config["vds_version"]}'
    )

    tmp_prefix_for_withdrawals: str = slugify(str(combiner_tmp_path))

    vds_paths: list[str] = []
    sg_ids_in_vds: list[str] = []
    new_sg_gvcfs: list[str] = []
    cohort_sgs: list[SequencingGroup] = cohort.get_sequencing_groups(only_active=True)
    cohort_sg_ids: list[str] = cohort.get_sequencing_group_ids(only_active=True)
    gvcf_external_header: str | None = None
    sequencing_group_names: list[str] | None = None

    for vds_id in config_retrieve(['combiner', 'vds_analysis_ids'], []):
        tmp_query_res, tmp_sg_ids_in_vds = get_vds_ids_output(vds_id)
        vds_paths.append(tmp_query_res)
        sg_ids_in_vds = sg_ids_in_vds + tmp_sg_ids_in_vds

    sgs_for_withdrawal: list[str] = [sg for sg in sg_ids_in_vds if sg not in cohort_sg_ids]

    if not config_retrieve(['combiner', 'merge_only_vds'], False):
        # Get SG IDs from the cohort object itself, rather than call Metamist.
        # Get VDS IDs first and filter out from this list
        new_sg_gvcfs = [str(sg.gvcf) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds]
        if new_sg_gvcfs:
            gvcf_external_header = new_sg_gvcfs[0]
            sequencing_group_names = [
                str(sg.id) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds
            ]

    if new_sg_gvcfs and len(new_sg_gvcfs) == 0 and len(vds_paths) <= 1 and not sgs_for_withdrawal:
        raise Exception

    combiner_job: PythonJob = _initalise_combiner_job(cohort=cohort)

    # Default to GRCh38 for reference if not specified
    combiner_job.call(
        _run,
        output_vds_path=output_vds_path,
        sequencing_type=workflow_config['sequencing_type'],
        tmp_prefix=f'{cohort.analysis_dataset.tmp_prefix()}combiner_temp_dir',
        genome_build=genome_build(),
        save_path=save_path,
        force_new_combiner=config_retrieve(['combiner', 'force_new_combiner']),
        sequencing_group_names=sequencing_group_names,
        gvcf_external_header=gvcf_external_header,
        gvcf_paths=new_sg_gvcfs,
        vds_paths=vds_paths,
        sgs_for_withdrawal=sgs_for_withdrawal,
        tmp_prefix_for_withdrawals=tmp_prefix_for_withdrawals,
        worker_memory=combiner_config.get('worker_memory', 'highmem'),
        worker_cores=combiner_config.get('worker_cores', 1),
        driver_memory=combiner_config.get('driver_memory', 'highmem'),
        driver_cores=combiner_config.get('driver_cores', 1),
    )

    return combiner_job


def _run(
    output_vds_path: str,
    sequencing_type: str,
    tmp_prefix: str,
    genome_build: str,
    save_path: str | None,
    sgs_for_withdrawal: list[str] | None = None,
    tmp_prefix_for_withdrawals: str | None = None,
    force_new_combiner: bool = False,
    sequencing_group_names: list[str] | None = None,
    gvcf_external_header: str | None = None,
    gvcf_paths: list[str] | None = None,
    vds_paths: list[str] | None = None,
    specific_intervals: list[str] | None = None,
    worker_memory: str = 'highmem',
    worker_cores: int = 1,
    driver_memory: str = 'highmem',
    driver_cores: int = 1,
) -> None:
    """
    Runs the combiner

    Scenarios, all predicated on output_vds_path not existing (otherwise this stage won't be scheduled) -
    1. no SGs to withdraw, combining to do
        - create a new combiner result, writing to the main output path
    2. SGs to withdraw, and combining to do
        - create a new combiner result, writing to the temporary output path
        - remove the withdrawn SGs from the temp VDS
    3. SGs to withdraw, no combining to do
        - don't run the combiner... this assumes vds_paths has only one element
        - filter samples from the result, then write to the main output path

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

    if vds_paths is None:
        vds_paths = []

    if sgs_for_withdrawal and tmp_prefix_for_withdrawals is None:
        raise ValueError('tmp_prefix_for_withdrawals must be set if sgs_for_withdrawal is set')

    import logging

    import hail as hl

    # set up a quick logger inside the job
    logging.basicConfig(level=logging.INFO)

    # init batch early, this method won't run if output_vds_path exists, so we have some work to do
    init_batch(
        worker_memory=worker_memory,
        driver_memory=driver_memory,
        driver_cores=driver_cores,
        worker_cores=worker_cores,
    )

    # do we need to do any combining?
    combining_to_do: bool = bool(gvcf_paths or len(vds_paths) > 1)

    # if we're going to do some combining, we need to set up the combiner instance
    if combining_to_do:
        if specific_intervals:
            logging.info(f'Using specific intervals: {specific_intervals}')
            intervals = hl.eval(
                [hl.parse_locus_interval(interval, reference_genome=genome_build) for interval in specific_intervals],
            )

        else:
            intervals = None

        # If we have SGs to withdraw, we write to a temp location, otherwise write straight out
        combiner_vds_output_path = tmp_prefix_for_withdrawals if sgs_for_withdrawal else output_vds_path

        # Load from save, if supplied
        if save_path:
            if force_new_combiner:
                logging.info(f'Combiner plan {save_path} will be ignored/written new')
            else:
                logging.info(f'Resuming combiner plan from {save_path}')

        combiner_instance: VariantDatasetCombiner = new_combiner(
            output_path=combiner_vds_output_path,
            save_path=save_path,
            gvcf_paths=gvcf_paths,
            gvcf_sample_names=sequencing_group_names,
            # Header must be used with gvcf_sample_names, otherwise gvcf_sample_names
            # will be ignored. The first gvcf path works fine as a header because it will
            # be only read until the last line that begins with "#":
            gvcf_external_header=gvcf_external_header,
            vds_paths=vds_paths,
            reference_genome=genome_build,
            temp_path=tmp_prefix,
            use_exome_default_intervals=sequencing_type == 'exome',
            use_genome_default_intervals=sequencing_type == 'genome',
            intervals=intervals,
            force=force_new_combiner,
            branch_factor=config_retrieve(
                ['combiner', 'branch_factor'],
                100,
            ),
            gvcf_batch_size=config_retrieve(
                ['combiner', 'gvcf_batch_size'],
                None,
            ),
        )

        combiner_instance.run()

    # do we need to run filter_samples?
    if sgs_for_withdrawal:
        # select which path to filter samples from - either the temp path if combining has occurred, or the only VDS
        path_to_read_from = tmp_prefix_for_withdrawals if combining_to_do else vds_paths[0]

        logging.info(f'There are {len(sgs_for_withdrawal)} sequencing groups to remove.')
        logging.info(f'Removing: {" ".join(sgs_for_withdrawal)}')
        filtered_vds: VariantDataset = hl.vds.filter_samples(
            vds=hl.vds.read_vds(path_to_read_from),
            samples=sgs_for_withdrawal,
            keep=False,
            remove_dead_alleles=True,
        )
        filtered_vds.write(output_vds_path)
