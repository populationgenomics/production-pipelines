import json
import logging
from typing import TYPE_CHECKING, Any, Final, Tuple

from cpg_utils import Path, to_path
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

    import hail as hl
    from hailtop.batch.job import PythonJob, PythonResult


HAIL_QUERY: Final = 'hail query'


# TODO, update analysis_meta here to pull the gvcf.type, and store this in metamist.
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
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.memory(combiner_config.get('memory'))
        j.storage(combiner_config.get('storage'))

        # set this job to be non-spot (i.e. non-preemptible)
        # previous issues with preemptible VMs led to multiple simultaneous QOB groups processing the same data
        j.spot(config_retrieve(['combiner', 'preemptible_vms'], False))

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(output_paths['vds']),
            sequencing_type=workflow_config['sequencing_type'],
            tmp_prefix=tmp_prefix,
            genome_build=genome_build(),
            save_path=output_paths['combiner_plan'],
            force_new_combiner=config_retrieve(['combiner', 'force_new_combiner']),
            sequencing_group_names=[
                str(sg.id) for sg in cohort_sgs if sg.gvcf is not None and sg.id not in sg_ids_in_vds
            ],
            gvcf_external_header=new_sg_gvcfs[0],
            gvcf_paths=new_sg_gvcfs,
            vds_paths=vds_paths,
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

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            if workflow_config.get(key):
                init_batch_args[key] = workflow_config[key]

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
        if site_only_version := config_retrieve(['large_cohort', 'output_versions', 'makesiteonly'], default=None):
            site_only_version = slugify(site_only_version)

        site_only_version = site_only_version or get_workflow().output_version

        return {
            'vcf': cohort.analysis_dataset.prefix() / get_workflow().name / site_only_version / 'siteonly.vcf.bgz',
            'tbi': cohort.analysis_dataset.prefix() / get_workflow().name / site_only_version / 'siteonly.vcf.bgz.tbi',
            'ht': cohort.analysis_dataset.prefix() / get_workflow().name / site_only_version / 'siteonly.ht',
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
                str(self.expected_outputs(cohort)['ht']),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=MakeSiteOnlyVcf)
class Vqsr(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if vqsr_version := config_retrieve(['large_cohort', 'output_versions', 'vqsr'], default=None):
            vqsr_version = slugify(vqsr_version)

        vqsr_version = vqsr_version or get_workflow().output_version
        return {
            'vcf': cohort.analysis_dataset.prefix() / get_workflow().name / vqsr_version / 'siteonly.vqsr.vcf.gz',
            'tbi': cohort.analysis_dataset.prefix() / get_workflow().name / vqsr_version / 'siteonly.vqsr.vcf.gz.tbi',
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
        if load_vqsr_version := config_retrieve(['large_cohort', 'output_versions', 'loadvqsr'], default=None):
            load_vqsr_version = slugify(load_vqsr_version)

        load_vqsr_version = load_vqsr_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / get_workflow().name / load_vqsr_version / 'vqsr.ht'

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


@stage(required_stages=[Combiner, Relatedness, Ancestry, MakeSiteOnlyVcf, LoadVqsr])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if frequencies_version := config_retrieve(['large_cohort', 'output_versions', 'frequencies'], default=None):
            frequencies_version = slugify(frequencies_version)

        frequencies_version = frequencies_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / get_workflow().name / frequencies_version / 'frequencies.ht'

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
                str(inputs.as_path(cohort, Ancestry, key='sample_qc_ht')),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(inputs.as_path(cohort, MakeSiteOnlyVcf, key='ht')),
                str(inputs.as_path(cohort, LoadVqsr)),
                str(self.expected_outputs(cohort)),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[Combiner])
class ShardVds(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if sharded_vds_version := config_retrieve(['large_cohort', 'output_versions', 'sharded_vds'], default=None):
            sharded_vds_version = slugify(sharded_vds_version)

        sharded_vds_version = sharded_vds_version or get_workflow().output_version
        return {
            f'{contig}': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / 'sharded_vds'
            / sharded_vds_version
            / f'{contig}.vds'
            for contig in [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import generate_coverage_table

        j = get_batch().new_job(
            'ShardVds',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                generate_coverage_table,
                generate_coverage_table.shard_vds.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                {k: str(v) for k, v in self.expected_outputs(cohort).items()},
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage()  # maybe not required?
class GenerateReferenceCoverageTable(CohortStage):
    """
    The `reference_ht` is a Table that contains a row for each locus coverage that should be
    computed on. It needs to be keyed by `locus`. The `reference_ht` can e.g. be
    created using `get_reference_ht`.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Reference coverage tables are created for each base in the region sequenced. This is an expensive operation
        and need only be done once per sequencing type (e.g. exome, genome) and exome capture method.
        Therefore, outputs are stored as a CPG-wide resource.
        """
        if config_retrieve(['workflow', 'sequencing_type']) == 'exome':
            ref_cov_version = config_retrieve(
                ['large_cohort', 'output_versions', 'exome_reference_coverage'],
                default=None,
            )
        else:
            ref_cov_version = config_retrieve(
                ['large_cohort', 'output_versions', 'genome_reference_coverage'],
                default=None,
            )
        contig_lengths_file = config_retrieve(['large_cohort', 'references', 'contig_lengths'], default=None)
        shard_size = config_retrieve(['large_cohort', 'interval_size'], default=500_000)

        with open(to_path(contig_lengths_file)) as f:
            contig_lengths: dict[str, int] = json.load(f)

        contigs = contig_lengths.keys()

        ref_cov_version = ref_cov_version or get_workflow().output_version
        return {
            f'{contig}_{start}_{end}': cohort.analysis_dataset.prefix()  # TODO: pick bucket to store in that's not cohort-specific
            / get_workflow().name
            / ref_cov_version
            / 'reference_coverage'
            / f'{contig}_reference_coverage_{start}_{end}.ht'
            for contig in contigs
            for start in range(1, contig_lengths[contig], shard_size)
            for end in [min(start + shard_size, contig_lengths[contig])]
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import generate_coverage_table

        jobs = []
        shard_size = config_retrieve(['large_cohort', 'interval_size'], default=500_000)
        congtig_lengths_file = config_retrieve(['large_cohort', 'references', 'contig_lengths'], default=None)

        outputs = self.expected_outputs(cohort)
        # TODO: detect when end of chromosome is reached and need to pass includes_end=True
        for shard, out_path in outputs.items():
            chrom, start, end = shard.split('_')
            j = get_batch().new_python_job(
                f'GenerateReferenceTable_{shard}',
                (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
            )
            j.image(image_path('cpg_workflows'))
            j.call(
                generate_coverage_table.generate_reference_coverage_ht,
                ref='GRCh38',
                chrom=chrom,
                start=int(start),
                end=int(end),
                shard_size=shard_size,
                out_path=str(out_path),
            )
            jobs.append(j)

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=[Combiner])
class GenerateCoverageTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if coverage_version := config_retrieve(['large_cohort', 'output_versions', 'coverage'], default=None):
            coverage_version = slugify(coverage_version)

        contig_lengths_file = config_retrieve(['large_cohort', 'references', 'contig_lengths'], default=None)
        shard_size = config_retrieve(['large_cohort', 'interval_size'], default=500_000)

        with open(to_path(contig_lengths_file)) as f:
            contig_lengths: dict[str, int] = json.load(f)

        contigs = contig_lengths.keys()

        coverage_version = coverage_version or get_workflow().output_version
        return {
            f'{contig}_{start}_{end}': cohort.analysis_dataset.prefix()  # TODO: pick bucket to store in that's not cohort-specific
            / get_workflow().name
            / coverage_version
            / 'split_coverage'
            / f'{contig}_coverage_{start}_{end}.ht'
            for contig in contigs
            for start in range(1, contig_lengths[contig], shard_size)
            for end in [min(start + shard_size, contig_lengths[contig])]
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import generate_coverage_table
        from cpg_workflows.large_cohort.generate_coverage_table import generate_intervals

        converage_jobs = []

        interval_size = config_retrieve(['large_cohort', 'interval_size'], default=500_000)

        coverage_table_paths = self.expected_outputs(cohort)

        for shard, coverage_output_path in coverage_table_paths.items():
            chrom, start, end = shard.split('_')
            j = get_batch().new_job(
                f'GenerateCoverageTable_{shard}',
                (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
            )
            j.image(config_retrieve(['workflow', 'driver_image']))
            j.command(
                query_command(
                    generate_coverage_table,
                    generate_coverage_table.run.__name__,
                    str(inputs.as_path(cohort, Combiner, key='vds')),
                    chrom,
                    int(start),
                    int(end),
                    str(coverage_output_path),
                    setup_gcp=True,
                ),
            )
            converage_jobs.append(j)

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=converage_jobs)


@stage(required_stages=[GenerateCoverageTable])
class MergeCoverageTables(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if coverage_version := config_retrieve(['large_cohort', 'output_versions', 'coverage'], default=None):
            coverage_version = slugify(coverage_version)

        coverage_version = coverage_version or get_workflow().output_version
        return (
            cohort.analysis_dataset.prefix()
            / get_workflow().name
            / coverage_version
            / 'merged_coverage'
            / 'merged_coverage.ht'
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import generate_coverage_table

        j = get_batch().new_job(
            'MergeCoverageTables',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                generate_coverage_table,
                generate_coverage_table.merge_coverage_tables.__name__,
                [str(v) for v in inputs.as_dict(cohort, GenerateCoverageTable).values()],
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


# @stage(required_stages=[Frequencies])
@stage()
class PrepareBrowserTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if browser_version := config_retrieve(['large_cohort', 'output_versions', 'preparebrowsertable'], default=None):
            browser_version = slugify(browser_version)

        browser_version = browser_version or get_workflow().output_version
        return {
            'browser': cohort.analysis_dataset.prefix() / get_workflow().name / browser_version / 'browser.ht',
            'exome_variants': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / browser_version
            / 'exome_variants.ht',
            'genome_variants': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / browser_version
            / 'genome_variants.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import browser_prepare

        j = get_batch().new_job(
            'PrepareBrowserTable',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        exome_freq_ht_path = config_retrieve(['large_cohort', 'output_versions', 'frequencies_exome'], default=None)
        genome_freq_ht_path = config_retrieve(['large_cohort', 'output_versions', 'frequencies_genome'], default=None)

        j.command(
            query_command(
                browser_prepare,
                browser_prepare.prepare_v4_variants.__name__,
                # hard-coding Frequencies tables for now
                exome_freq_ht_path,
                genome_freq_ht_path,
                str(self.expected_outputs(cohort)['browser']),
                str(self.expected_outputs(cohort)['exome_variants']),
                str(self.expected_outputs(cohort)['genome_variants']),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
