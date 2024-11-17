from cpg_utils import Path
from cpg_utils.config import config_retrieve, genome_build, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.resources import joint_calling_scatter_count
from cpg_workflows.targets import Cohort, MultiCohort, SequencingGroup
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    CohortStage,
    MultiCohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
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


@stage(analysis_type='combiner')
class GVCFCombiner(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.prefix / f'{multicohort.name}.vds'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import combiner

        output_vds_path: Path = self.expected_outputs(multicohort)

        # create these as empty lists instead of None, they have the same truthiness
        vds_path: str | None = None
        sg_ids_in_vds: list[str] = []

        if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
            vds_path = existing_vds_analysis_entry['output']
            sg_ids_in_vds = [sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']]

        new_sg_gvcfs: list[str] = [
            str(sg.gvcf)
            for sg in multicohort.get_sequencing_groups()
            if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        ]

        if not new_sg_gvcfs:
            return self.make_outputs(multicohort, output_vds_path)

        j = get_batch().new_python_job('Combiner', self.get_job_attrs())
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.memory(config_retrieve(['workflow', 'memory']))
        j.storage(config_retrieve(['workflow', 'storage']))

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(output_vds_path),
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_paths=[vds_path],
        )

        return self.make_outputs(multicohort, output_vds_path, j)


@stage(required_stages=[GVCFCombiner], analysis_type='matrixtable', analysis_keys=['mt', 'vcf_dir'])
class DenseMTFromVDS(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        return {
            'mt': self.prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF
            'vcf_dir': str(self.prefix / f'{multicohort.name}.vcf.bgz'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        output = self.expected_outputs(multicohort)

        j = get_batch().new_job('Dense Subset')
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.command(
            'mt_from_vds '
            f'--input {str(inputs.as_path(multicohort, GVCFCombiner))} '
            f'--output {str(output["mt"])} '
            f'--sites_only {output["vcf_dir"]}'
            f'--partitions {joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))}',
        )
        return self.make_outputs(multicohort, output, [j])


@stage(required_stages=DenseMTFromVDS)
class Vqsr(CohortStage):
    """
    Exciting. A chance to try this out.

    The parallel export of VCFs exports a collection of VCFs and indexes,
    and a manifest file.
    By reading this manifest file as a python Job, we can find the single VCF by name relevant to each interval
    We can localise it in a second job, pass it into the next job to run VQSR on it.
    We can use the same logic to grab the same VCFs for VEP
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'siteonly.vqsr.vcf.gz',
            'tbi': self.tmp_prefix / 'siteonly.vqsr.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.jobs import vqsr

        vcf_path = inputs.as_path(cohort, MakeSiteOnlyVcf, key='vcf')
        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=vcf_path,
            gvcf_count=len(cohort.get_sequencing_groups()),
            out_path=self.expected_outputs(cohort)['vcf'],
            tmp_prefix=self.tmp_prefix,
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=Vqsr)
class LoadVqsr(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return get_workflow().prefix / 'vqsr.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import load_vqsr

        j = get_batch().new_job('LoadVqsr', self.get_job_attrs())
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                load_vqsr,
                load_vqsr.run.__name__,
                str(inputs.as_path(cohort, Vqsr, key='vcf')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
