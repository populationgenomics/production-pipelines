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
    workflow_config = config_retrieve('workflow')
    combiner_config = config_retrieve('combiner')

    tmp_prefix_for_withdrawals: str = slugify(
        f'{cohort.analysis_dataset.tmp_prefix}/{workflow_config["cohort"]}-{workflow_config["sequencing_type"]}-{combiner_config["vds_version"]}',
    )

    vds_paths: list[str] = []
    sg_ids_in_vds: list[str] = []
    sgs_for_withdrawal: list[str] = []
    new_sg_gvcfs: list[str] = []
    cohort_sgs: list[SequencingGroup] = cohort.get_sequencing_groups(only_active=True)
    cohort_sg_ids: list[str] = cohort.get_sequencing_group_ids(only_active=True)

    if combiner_config.get('vds_analysis_ids', None) is not None:
        for vds_id in combiner_config['vds_analysis_ids']:
            tmp_query_res, tmp_sg_ids_in_vds = get_vds_ids_output(vds_id)
            vds_paths.append(tmp_query_res)
            sg_ids_in_vds = sg_ids_in_vds + tmp_sg_ids_in_vds
            sgs_for_withdrawal = [sg for sg in sg_ids_in_vds if sg not in cohort_sg_ids]

    if combiner_config.get('merge_only_vds', False) is not True:
        # Get SG IDs from the cohort object itself, rather than call Metamist.
        # Get VDS IDs first and filter out from this list
        new_sg_gvcfs = [str(sg.gvcf) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds]

    if new_sg_gvcfs and len(new_sg_gvcfs) == 0 and len(vds_paths) <= 1:
        raise Exception

    sequencing_group_names: list[str] = [
        str(sg.id) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds
    ]

    combiner_job: PythonJob = _initalise_combiner_job(cohort=cohort)

    # Default to GRCh38 for reference if not specified
    combiner_job.call(
        _run,
        output_vds_path=output_vds_path,
        sequencing_type=workflow_config['sequencing_type'],
        tmp_prefix=f'{cohort.analysis_dataset.tmp_prefix}/combiner_temp_dir',
        genome_build=genome_build(),
        save_path=save_path,
        force_new_combiner=config_retrieve(['combiner', 'force_new_combiner']),
        sequencing_group_names=sequencing_group_names,
        gvcf_external_header=new_sg_gvcfs[0],
        gvcf_paths=new_sg_gvcfs,
        vds_paths=vds_paths,
        sgs_for_withdrawal=sgs_for_withdrawal,
        tmp_prefix_for_withdrawals=tmp_prefix_for_withdrawals,
    )

    return combiner_job


def _run(
    output_vds_path: str,
    sequencing_type: str,
    tmp_prefix: str,
    genome_build: str,
    save_path: str | None,
    sgs_for_withdrawl: list[str] | None,
    tmp_prefix_for_withdrawals: str,
    force_new_combiner: bool = False,
    sequencing_group_names: list[str] | None = None,
    gvcf_external_header: str | None = None,
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

    if not can_reuse(to_path(output_vds_path)) or sgs_for_withdrawl:
        init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)

        if specific_intervals:
            logging.info(f'Using specific intervals: {specific_intervals}')

            intervals = hl.eval(
                [hl.parse_locus_interval(interval, reference_genome=genome_build) for interval in specific_intervals],
            )

        else:
            intervals = None

        if not sgs_for_withdrawl:
            # Load from save, if supplied
            if save_path:
                if force_new_combiner:
                    logging.info(f'Combiner plan {save_path} will be ignored/written new')
                else:
                    logging.info(f'Resuming combiner plan from {save_path}')

            combiner: VariantDatasetCombiner = new_combiner(
                output_path=output_vds_path,
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

            combiner.run()
        else:
            combiner = new_combiner(
                output_path=tmp_prefix_for_withdrawals,
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

            combiner.run()

            filtered_vds: VariantDataset = hl.vds.filter_samples(
                vds=tmp_prefix_for_withdrawals,
                samples=sgs_for_withdrawl,
                keep=False,
                remove_dead_alleles=True,
            )
            filtered_vds.write(output_vds_path)
