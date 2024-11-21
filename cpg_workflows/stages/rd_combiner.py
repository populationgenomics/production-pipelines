from cpg_utils import Path
from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import MultiCohort
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    MultiCohortStage,
    StageInput,
    StageOutput,
    stage,
)
from metamist.graphql import gql, query

LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                meta
                output
                sequencingGroups {
                    id
                }
                timestampCompleted
            }
        }
    }
""",
)


def query_for_latest_vds(dataset: str, entry_type: str = 'combiner') -> dict | None:
    """
    query for the latest analysis of type entry_type for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest analysis
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': entry_type})

    analyses_by_date = {}

    for analysis in result['project']['analyses']:
        if analysis['output'] and (
            analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type'])
        ):
            analyses_by_date[analysis['timestampCompleted']] = analysis

    if not analyses_by_date:
        get_logger(__file__).warning(f'No analysis of type {entry_type} found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort as strings
    return analyses_by_date[sorted(analyses_by_date)[-1]]


@stage(analysis_type='combiner', analysis_keys=['vds'])
class GVCFCombiner(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'combiner_plan': str(self.tmp_prefix / 'combiner_plan.json'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import combiner

        outputs: dict[str, str | Path] = self.expected_outputs(multicohort)

        # create these as empty lists instead of None, they have the same truthiness
        vds_path: str | None = None
        sg_ids_in_vds: list[str] = []

        if config_retrieve(['workflow', 'check_for_existing_vds'], True):
            # check for existing VDS
            get_logger(__file__).info('Checking for existing VDS')
            if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
                vds_path = existing_vds_analysis_entry['output']
                sg_ids_in_vds = [sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']]

        new_sg_gvcfs: list[str] = [
            str(sg.gvcf)
            for sg in multicohort.get_sequencing_groups()
            if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        ]

        if not new_sg_gvcfs:
            get_logger(__file__).info('No GVCFs to combine')
            get_logger(__file__).info(f'Checking if VDS exists: {outputs["vds"]}: {outputs["vds"].exists()}')  # type: ignore
            return self.make_outputs(multicohort, outputs)

        j = get_batch().new_python_job('Combiner', self.get_job_attrs())
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.memory(config_retrieve(['combiner', 'memory']))
        j.storage(config_retrieve(['combiner', 'storage']))

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(outputs['vds']),
            save_path=outputs['combiner_plan'],
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_paths=[vds_path] if vds_path else None,
        )

        return self.make_outputs(multicohort, outputs, j)


@stage(required_stages=[GVCFCombiner], analysis_type='matrixtable', analysis_keys=['mt'])
class DenseMTFromVDS(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        return {
            'mt': self.prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF
            'vcf_dir': str(self.prefix / f'{multicohort.name}.vcf.bgz'),
            # this will be the file which contains the name of all fragments
            'shard_manifest': str(self.prefix / f'{multicohort.name}.vcf.bgz' / 'shard_manifest.txt'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        output = self.expected_outputs(multicohort)

        j = get_batch().new_job('Dense Subset')
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.command(
            'mt_from_vds '
            f'--input {str(inputs.as_dict(multicohort, GVCFCombiner)["vds"])} '
            f'--output {str(output["mt"])} '
            f'--sites_only {output["vcf_dir"]}',
        )
        return self.make_outputs(multicohort, output, [j])
