import logging
from functools import cache

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import ExpectedResultT, slugify, tshirt_mt_sizing
from cpg_workflows.workflow import CohortStage, StageInput, StageOutput, get_workflow, stage

from .genotype import Genotype


@cache
def relatedness_version() -> str:
    """
    generate the relatedness version
    """
    if not (relatedness_str := config_retrieve(['large_cohort', 'output_versions', 'relatedness'], default=None)):
        return get_workflow().output_version
    return slugify(relatedness_str)


@cache
def sample_qc_version() -> str:
    """
    generate the sample_qc version
    """
    if not (sample_qc_str := config_retrieve(['large_cohort', 'output_versions', 'sample_qc'], default=None)):
        return get_workflow().output_version
    return slugify(sample_qc_str)


@cache
def vds_version() -> str:
    if not (vds_version_str := get_config()['workflow'].get('vds_version')):
        return get_workflow().output_version
    vds_version_str = slugify(vds_version_str)
    if not vds_version_str.startswith('v'):
        vds_version_str = f'v{vds_version_str}'
    return vds_version_str


@cache
def dense_subset_version() -> str:
    """
    generate the dense_subset version
    """
    if not (dense_subset_str := config_retrieve(['large_cohort', 'output_versions', 'dense_subset'], default=None)):
        return get_workflow().output_version
    return slugify(dense_subset_str)


@stage(required_stages=[Genotype])
class Combiner(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return cohort.analysis_dataset.prefix() / 'vds' / f'{vds_version()}.vds'

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
        return cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version() / 'sample_qc.ht'

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
        return cohort.analysis_dataset.prefix() / get_workflow().name / dense_subset_version() / 'dense_subset.mt'

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


@stage(required_stages=DenseSubset, analysis_keys=['relatedness_ht'], analysis_type='custom')
class RelatednessPCRelate(CohortStage):
    """
    This is the first step of the Relatedness Stage, without Dataproc
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'relatedness_ht': (
                cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version() / 'relatedness.ht'
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(cohort)
        dense_mt_path = str(inputs.as_path(cohort, DenseSubset))
        dense_mt_name = dense_mt_path.split('/')[-1]

        # estimate required storage - we create a new HT, so provision double the storage
        t_shirt_size_value = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        ).split('Gi')[0]
        required_storage_value = int(t_shirt_size_value) * 2
        required_storage = f'{required_storage_value}Gi'

        # create a job
        job = get_batch().new_job(f'Run Relatedness PCRelate stage: {cohort.name}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command('set -eux pipefail')

        # localise the Dense MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {dense_mt_path} $BATCH_TMPDIR')

        # how many cores do we need?
        logging.info(f'Required storage: {required_storage}')
        required_cpu: int = config_retrieve(['RelatednessPCRelate', 'cores'], 8)
        job.cpu(required_cpu).storage(required_storage).memory('highmem')

        job.command(
            f'relatedness_pcrelate '
            f'--dense_mt "${{BATCH_TMPDIR}}/{dense_mt_name}" '
            f'--out {job.output} '
            f'--checkpoint "${{BATCH_TMPDIR}}"',
        )

        # delocalise the output HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.output} {str(outputs["relatedness_ht"])}')

        return self.make_outputs(cohort, outputs, job)


@stage(
    required_stages=[SampleQC, RelatednessPCRelate],
    analysis_keys=['relateds_to_drop'],
    analysis_type='custom',
)
class RelatednessFlag(CohortStage):

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'relateds_to_drop': (
                cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version() / 'relateds_to_drop.ht'
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(cohort)

        # estimate required storage - I've got nothing to base this on...
        t_shirt_size_value = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        ).split('Gi')[0]
        required_storage_value = int(t_shirt_size_value) * 2
        required_storage = f'{required_storage_value}Gi'

        # create a job
        job = get_batch().new_job(f'Run Relatedness Sample Flagging stage: {cohort.name}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command('set -eux pipefail')

        # localise the Relatedness HT
        prcrelate_mt_path = str(inputs.as_path(cohort, RelatednessPCRelate, 'relatedness_ht'))
        prcrelate_mt_name = prcrelate_mt_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {prcrelate_mt_path} $BATCH_TMPDIR')

        # localise the Sample QC HT MT
        sample_qc_ht_path = str(inputs.as_path(cohort, SampleQC))
        sample_qc_ht_name = sample_qc_ht_path.split('/')[-1]

        # localise the Dense MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {sample_qc_ht_path} $BATCH_TMPDIR')

        # how many cores do we need?
        required_cpu: int = config_retrieve(['RelatednessFlag', 'cores'], 8)
        job.cpu(required_cpu).storage(required_storage).memory('highmem')

        job.command(
            'relatedness_flag '
            f'--relatedness "${{BATCH_TMPDIR}}/{prcrelate_mt_name}" '
            f'--qc "${{BATCH_TMPDIR}}/{sample_qc_ht_name}" '
            f'--out {job.output} '
            f'--checkpoint "${{BATCH_TMPDIR}}" ',
        )

        # delocalise the output HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.output} {str(outputs["relateds_to_drop"])}')

        return self.make_outputs(cohort, outputs, job)


@stage(required_stages=[DenseSubset, SampleQC])
class AncestryAddBackground(CohortStage):
    """
    optional stage, runs if we need to annotate with background datasets
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'bg_qc_ht': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / sample_qc_version()
            / 'bg_sample_qc.ht',
            'bg_dense_mt': cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version() / 'bg_dense.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        # no datasets, no running
        if not config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            return self.make_outputs(target=cohort)

        # inputs
        dense_mt_path = inputs.as_path(cohort, DenseSubset)
        sample_qc_ht_path = inputs.as_path(cohort, SampleQC)

        # expected outputs
        outputs = self.expected_outputs(cohort)

        # job runs QOB, shouldn't need many resources itself
        job = get_batch().new_job('Add Ancestry Background')
        job.storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))
        job.command(
            f'ancestry_add_background '
            f'--qc_in {str(sample_qc_ht_path)} '
            f'--dense_in {str(dense_mt_path)} '
            f'--qc_out {str(outputs["bg_qc_ht"])} '
            f'--dense_out {str(outputs["bg_dense_mt"])} '
            f'--tmp {str(self.tmp_prefix)} ',
        )

        return self.make_outputs(target=cohort, jobs=job, data=outputs)


@stage(required_stages=[SampleQC, RelatednessFlag, DenseSubset, AncestryAddBackground])
class AncestryPCA(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'scores': ancestry_prefix / 'scores.ht',
            'eigenvalues': ancestry_prefix / 'eigenvalues.ht',
            'loadings': ancestry_prefix / 'loadings.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        # decide which stage to pull input from
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            dense_mt_path = str(inputs.as_path(cohort, AncestryAddBackground, 'bg_dense_mt'))

        else:
            dense_mt_path = str(inputs.as_path(cohort, DenseSubset))
            # sample_qc_ht_path = inputs.as_path(cohort, SampleQC)

        # negligible amount of storage for this one
        related = str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop'))

        # big job
        job = get_batch().new_bash_job('Run Ancestry PCA')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('50gi').memory('highmem').cpu(config_retrieve(['Ancestry', 'cores'], 8))

        # localise the Dense MT
        dense_name = dense_mt_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {dense_mt_path} $BATCH_TMPDIR')

        # localise the relatedness HT
        related_name = related.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {related} $BATCH_TMPDIR')

        # run the command
        job.command(
            'ancestry_pca '
            f'--dense_mt "${{BATCH_TMPDIR}}/{dense_name}" '
            f'--related "${{BATCH_TMPDIR}}/{related_name}" '
            f'--scores_out {job.scores} '
            f'--eigen_out {job.eigen} '
            f'--loadings_out {job.loadings} ',
        )

        # copy 3 tables out
        outputs = self.expected_outputs(cohort)
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.scores} {str(outputs["scores"])}')
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.eigen} {str(outputs["eigenvalues"])}')
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.loadings} {str(outputs["loadings"])}')

        return self.make_outputs(target=cohort, data=outputs, jobs=job)


@stage(required_stages=[AncestryPCA, SampleQC, AncestryAddBackground])
class AncestryInfer(CohortStage):

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:

        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'inferred_pop': ancestry_prefix / 'inferred_pop.ht',
            'sample_qc_ht': ancestry_prefix / 'sample_qc_ht.ht',
            'pickled_model': ancestry_prefix / 'pop.RF_fit.pickle',
            'inference_txt': ancestry_prefix / 'RF_pop_assignments.txt.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            sample_qc_ht_path = str(inputs.as_path(cohort, AncestryAddBackground, 'bg_qc_ht'))
        else:
            sample_qc_ht_path = str(inputs.as_path(cohort, SampleQC))
        scores_table = str(inputs.as_path(cohort, AncestryPCA, 'scores'))

        # big job
        job = get_batch().new_bash_job('Run Ancestry PCA')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('50gi').memory('highmem').cpu(config_retrieve(['Ancestry', 'cores'], 8))

        sample_qc_name = sample_qc_ht_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {sample_qc_ht_path} $BATCH_TMPDIR')

        scores_name = scores_table.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {scores_table} $BATCH_TMPDIR')

        # run the command
        job.command(
            'ancestry_infer_labels '
            f'--scores_ht "${{BATCH_TMPDIR}}/{scores_name}" '
            f'--qc_in "${{BATCH_TMPDIR}}/{sample_qc_name}" '
            f'--qc_out {job.qc} '
            f'--ht_out {job.ht} '
            f'--pickle_out {job.pickle} '
            f'--txt_out {job.txt} ',
        )

        outputs = self.expected_outputs(cohort)

        # copy out the scores HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.ht} {str(outputs["inferred_pop"])}')

        # copy out the sample QC HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.qc} {str(outputs["sample_qc_ht"])}')

        # and copy the other single files
        get_batch().write_output(job.pickle, str(outputs['pickled_model']))
        get_batch().write_output(job.txt, str(outputs['inference_txt']))

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=[SampleQC, DenseSubset, RelatednessFlag])
class Ancestry(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'scores': ancestry_prefix / 'scores.ht',
            'eigenvalues': ancestry_prefix / 'eigenvalues.ht',
            'loadings': ancestry_prefix / 'loadings.ht',
            'inferred_pop': ancestry_prefix / 'inferred_pop.ht',
            'sample_qc_ht': ancestry_prefix / 'sample_qc_ht.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_pca import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop'),
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


@stage(required_stages=[SampleQC, AncestryInfer, RelatednessFlag])
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
                relateds_to_drop_ht_path=inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop'),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner, SampleQC, RelatednessFlag])
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
                str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop')),
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
    def expected_outputs(self, cohort: Cohort) -> Path:
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


@stage(required_stages=[Combiner, SampleQC, RelatednessFlag])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
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
                str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
