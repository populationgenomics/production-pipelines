import logging
from typing import TYPE_CHECKING, Any, Final, Tuple

from cpg_utils import Path
from cpg_utils.config import config_retrieve, genome_build, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.utils import slugify
from cpg_workflows.workflow import (
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

if TYPE_CHECKING:
    from graphql import DocumentNode

    from hailtop.batch.job import PythonJob


HAIL_QUERY: Final = 'hail query'


@stage(analysis_type='combiner', analysis_keys=['vds'])
class Combiner(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Any]:
        workflow_config = config_retrieve('workflow')
        combiner_config = config_retrieve('combiner')
        output_vds_name: str = slugify(
            f"{cohort.name}-{workflow_config['sequencing_type']}-{combiner_config['vds_version']}",
        )

        # include the list of all VDS IDs in the plan name
        if vds_ids := config_retrieve(['combiner', 'vds_analysis_ids']):
            ids_list_as_string: str = '_'.join([str(id) for id in sorted(vds_ids)])
            combiner_plan_name: str = f'combiner_{ids_list_as_string}'
        else:
            combiner_plan_name = f'combiner-{cohort.name}'
        return {
            'vds': cohort.analysis_dataset.prefix() / 'vds' / f'{output_vds_name}.vds',
            'combiner_plan': str(self.get_stage_cohort_prefix(cohort, 'tmp') / f'{combiner_plan_name}.json'),
        }

    def get_vds_ids_output(self, vds_id: int) -> Tuple[str, list[str]]:
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

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        # Can't import it before all configs are set:
        from cpg_workflows.large_cohort import combiner

        workflow_config = config_retrieve('workflow')
        combiner_config = config_retrieve('combiner')

        output_paths = self.expected_outputs(cohort)
        tmp_prefix = slugify(
            f"{self.tmp_prefix}/{workflow_config['cohort']}-{workflow_config['sequencing_type']}-{combiner_config['vds_version']}",
        )

        # create these as empty lists instead of None, they have the same truthiness
        vds_paths: list[str] = []
        sg_ids_in_vds: list[str] = []
        new_sg_gvcfs: list[str] = []

        if combiner_config.get('vds_analysis_ids', None) is not None:
            for vds_id in combiner_config['vds_analysis_ids']:
                tmp_query_res, tmp_sg_ids_in_vds = self.get_vds_ids_output(vds_id)
                vds_paths.append(tmp_query_res)
                sg_ids_in_vds = sg_ids_in_vds + tmp_sg_ids_in_vds

        if combiner_config.get('merge_only_vds', False) is not True:
            # Get SG IDs from the cohort object itself, rather than call Metamist.
            # Get VDS IDs first and filter out from this list
            cohort_sgs: list[SequencingGroup] = cohort.get_sequencing_groups(only_active=True)
            new_sg_gvcfs = [str(sg.gvcf) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds]

        if new_sg_gvcfs and len(new_sg_gvcfs) == 0 and len(vds_paths) <= 1:
            return self.make_outputs(cohort, output_paths)

        j: PythonJob = get_batch().new_python_job('Combiner', (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY})
        j.image(image_path('cpg_workflows'))
        j.memory(combiner_config.get('memory'))
        j.storage(combiner_config.get('storage'))

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(output_paths['vds']),
            sequencing_type=workflow_config['sequencing_type'],
            tmp_prefix=tmp_prefix,
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_paths=vds_paths,
            save_path=output_paths['combiner_plan'],
            force_new_combiner=config_retrieve(['combiner', 'force_new_combiner']),
        )

        return self.make_outputs(cohort, output_paths, [j])


@stage(required_stages=[Combiner])
class SampleQC(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if sample_qc_version := config_retrieve(['large_cohort', 'output_versions', 'sample_qc'], default=None):
            sample_qc_version = slugify(sample_qc_version)

        sample_qc_version = sample_qc_version or get_workflow().output_version
        sample_qc_path = cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version / 'sample_qc.ht'
        return sample_qc_path

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import sample_qc

        j = get_batch().new_job(
            'Sample QC',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        if 'driver_cores' in workflow_config:
            init_batch_args['driver_cores'] = workflow_config['driver_cores']

        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                sample_qc,
                sample_qc.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(self.expected_outputs(cohort)),
                str(self.tmp_prefix),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class DenseSubset(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if dense_subset_version := config_retrieve(['large_cohort', 'output_versions', 'dense_subset'], default=None):
            dense_subset_version = slugify(dense_subset_version)

        dense_subset_version = dense_subset_version or get_workflow().output_version
        dense_subset_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / dense_subset_version / 'dense_subset.mt'
        )
        return dense_subset_path

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import dense_subset

        j = get_batch().new_job(
            'Dense Subset',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                dense_subset,
                dense_subset.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, DenseSubset])
class Relatedness(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if relatedness_version := config_retrieve(['large_cohort', 'output_versions', 'relatedness'], default=None):
            relatedness_version = slugify(relatedness_version)

        relatedness_version = relatedness_version or get_workflow().output_version
        relatedness_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version / 'relatedness.ht'
        )
        relatedness_to_drop_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version / 'relateds_to_drop.ht'
        )

        return dict(
            relatedness=relatedness_path,
            relateds_to_drop=relatedness_to_drop_path,
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.relatedness import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                out_relatedness_ht_path=self.expected_outputs(cohort)['relatedness'],
                out_relateds_to_drop_ht_path=self.expected_outputs(cohort)['relateds_to_drop'],
                tmp_prefix=self.tmp_prefix,
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, DenseSubset, Relatedness])
class Ancestry(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return dict(
            scores=get_workflow().prefix / 'ancestry' / 'scores.ht',
            eigenvalues=get_workflow().prefix / 'ancestry' / 'eigenvalues.ht',
            loadings=get_workflow().prefix / 'ancestry' / 'loadings.ht',
            inferred_pop=get_workflow().prefix / 'ancestry' / 'inferred_pop.ht',
            sample_qc_ht=get_workflow().prefix / 'ancestry' / 'sample_qc_ht.ht',
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_pca import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort,
                    Relatedness,
                    key='relateds_to_drop',
                ),
                tmp_prefix=self.tmp_prefix,
                out_scores_ht_path=self.expected_outputs(cohort)['scores'],
                out_eigenvalues_ht_path=self.expected_outputs(cohort)['eigenvalues'],
                out_loadings_ht_path=self.expected_outputs(cohort)['loadings'],
                out_inferred_pop_ht_path=self.expected_outputs(cohort)['inferred_pop'],
                out_sample_qc_ht_path=self.expected_outputs(cohort)['sample_qc_ht'],
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, Ancestry, Relatedness])
class AncestryPlots(CohortStage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.out_prefix = get_workflow().web_prefix / 'ancestry'
        self.out_fname_pattern = '{scope}_pc{pci}_{pca_suffix}.{ext}'

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        n_pcs = get_config()['large_cohort']['n_pcs']
        # if there is a pca_plot_name given, add this to the output name
        plot_name = get_config()['large_cohort'].get('pca_plot_name')
        pca_suffix = ''
        if plot_name:
            pca_suffix = plot_name.replace('-', '_')
        return {
            str(pc_num): self.out_prefix
            / self.out_fname_pattern.format(
                scope='dataset',
                pci=pc_num,
                pca_suffix=pca_suffix,
                ext='html',
            )
            for pc_num in range(1, n_pcs)
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_plots import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                out_path_pattern=self.out_prefix / self.out_fname_pattern,
                sample_qc_ht_path=inputs.as_path(cohort, Ancestry, key='sample_qc_ht'),
                scores_ht_path=inputs.as_path(cohort, Ancestry, key='scores'),
                eigenvalues_ht_path=inputs.as_path(cohort, Ancestry, key='eigenvalues'),
                loadings_ht_path=inputs.as_path(cohort, Ancestry, key='loadings'),
                inferred_pop_ht_path=inputs.as_path(
                    cohort,
                    Ancestry,
                    key='inferred_pop',
                ),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort,
                    Relatedness,
                    key='relateds_to_drop',
                ),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner, SampleQC, Relatedness])
class MakeSiteOnlyVcf(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'siteonly.vcf.bgz',
            'tbi': self.tmp_prefix / 'siteonly.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import site_only_vcf

        j = get_batch().new_job(
            'MakeSiteOnlyVcf',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )

        sitesvcf_config = config_retrieve('sitesvcf')

        j.image(image_path('cpg_workflows'))
        j.memory(sitesvcf_config.get('memory', '4Gi'))
        j.storage(sitesvcf_config.get('storage', '5Gi'))

        init_batch_args: dict[str, str | int] = {}
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if sitesvcf_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        if 'driver_cores' in sitesvcf_config:
            init_batch_args['driver_cores'] = sitesvcf_config['driver_cores']

        j.command(
            query_command(
                site_only_vcf,
                site_only_vcf.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)['vcf']),
                str(self.tmp_prefix),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=MakeSiteOnlyVcf)
class Vqsr(CohortStage):
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

        j = get_batch().new_job(
            'LoadVqsr',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
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


@stage(required_stages=[Combiner, SampleQC, Relatedness])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return get_workflow().prefix / 'frequencies.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import frequencies

        j = get_batch().new_job(
            'Frequencies',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        if 'driver_cores' in workflow_config:
            init_batch_args['driver_cores'] = workflow_config['driver_cores']

        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                frequencies,
                frequencies.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
