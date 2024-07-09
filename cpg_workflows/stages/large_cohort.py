import logging

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import slugify
from cpg_workflows.workflow import (
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)

from .genotype import Genotype


@stage(required_stages=[Genotype])
class Combiner(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if vds_version := get_config()['workflow'].get('vds_version'):
            vds_version = slugify(vds_version)
            if not vds_version.startswith('v'):
                vds_version = f'v{vds_version}'

        vds_version = vds_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / 'vds' / f'{vds_version}.vds'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        # Can't import it before all configs are set:
        from cpg_workflows.large_cohort import combiner

        j = get_batch().new_job('Combiner', (self.get_job_attrs() or {}) | {'tool': 'hail query'})

        init_batch_args: dict[str, str | int] = {}
        config = config_retrieve('workflow')

        if config.get('highmem_workers'):
            init_batch_args['worker_memory'] = 'highmem'
        if config.get('highmem_drivers'):
            init_batch_args['driver_memory'] = 'highmem'
        if 'driver_cores' in config:
            init_batch_args['driver_cores'] = config['driver_cores']
        if not init_batch_args:
            logging.warning(
                "None of 'highmem_workers', 'highmem_drivers', or 'driver_cores' were specified in the config. If you're getting OOM errors, ensure these are included in the config.",
            )

        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                combiner,
                combiner.run.__name__,
                str(self.expected_outputs(cohort)),
                str(self.tmp_prefix),
                setup_gcp=True,
                init_batch_args=init_batch_args,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


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
            (self.get_job_attrs() or {}) | {'tool': 'hail query'},
        )
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                sample_qc,
                sample_qc.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(self.expected_outputs(cohort)),
                str(self.tmp_prefix),
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
            (self.get_job_attrs() or {}) | {'tool': 'hail query'},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                dense_subset,
                dense_subset.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
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
            (self.get_job_attrs() or {}) | {'tool': 'hail query'},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                site_only_vcf,
                site_only_vcf.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)['vcf']),
                str(self.tmp_prefix),
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
            (self.get_job_attrs() or {}) | {'tool': 'hail query'},
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
            (self.get_job_attrs() or {}) | {'tool': 'hail query'},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                frequencies,
                frequencies.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
