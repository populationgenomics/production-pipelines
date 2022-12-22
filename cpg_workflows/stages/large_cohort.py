from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import slugify
from cpg_workflows.workflow import (
    stage,
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
)
from .genotype import Genotype
from .. import get_batch


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
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.combiner import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                out_vds_path=self.expected_outputs(cohort),
                tmp_prefix=self.tmp_prefix,
            ),
            autoscaling_policy=(
                get_config()['hail']
                .get('dataproc', {})
                .get('combiner_autoscaling_policy')
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class SampleQC(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return get_workflow().prefix / 'sample_qc.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.sample_qc import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                vds_path=inputs.as_path(cohort, Combiner),
                out_sample_qc_ht_path=self.expected_outputs(cohort),
                tmp_prefix=self.tmp_prefix,
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class DenseSubset(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return get_workflow().prefix / 'dense_subset.mt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.dense_subset import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                vds_path=inputs.as_path(cohort, Combiner),
                out_dense_mt_path=self.expected_outputs(cohort),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, DenseSubset])
class Relatedness(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return dict(
            relatedness=get_workflow().prefix / 'relatedness.ht',
            relateds_to_drop=get_workflow().prefix / 'relateds_to_drop.ht',
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
                out_relateds_to_drop_ht_path=self.expected_outputs(cohort)[
                    'relateds_to_drop'
                ],
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
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.ancestry_pca import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort, Relatedness, key='relateds_to_drop'
                ),
                tmp_prefix=self.tmp_prefix,
                out_scores_ht_path=self.expected_outputs(cohort)['scores'],
                out_eigenvalues_ht_path=self.expected_outputs(cohort)['eigenvalues'],
                out_loadings_ht_path=self.expected_outputs(cohort)['loadings'],
                out_inferred_pop_ht_path=self.expected_outputs(cohort)['inferred_pop'],
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, Ancestry])
class AncestryPlots(CohortStage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.out_prefix = get_workflow().web_prefix / 'ancestry'
        self.out_fname_pattern = '{scope}_pc{pci}.{ext}'

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.out_prefix / self.out_fname_pattern.format(
            scope='dataset', pci=1, ext='html'
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.ancestry_plots import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                out_path_pattern=self.out_prefix / self.out_fname_pattern,
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                scores_ht_path=inputs.as_path(cohort, Ancestry, key='scores'),
                eigenvalues_ht_path=inputs.as_path(cohort, Ancestry, key='eigenvalues'),
                loadings_ht_path=inputs.as_path(cohort, Ancestry, key='loadings'),
                inferred_pop_ht_path=inputs.as_path(
                    cohort, Ancestry, key='inferred_pop'
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
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.site_only_vcf import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                vds_path=inputs.as_path(cohort, Combiner),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort, Relatedness, key='relateds_to_drop'
                ),
                out_vcf_path=self.expected_outputs(cohort)['vcf'],
                tmp_prefix=self.tmp_prefix,
            ),
            depends_on=inputs.get_jobs(cohort),
            # hl.export_vcf() uses non-preemptible workers' disk to merge VCF files.
            # 10 samples take 2.3G, 400 samples take 60G, which roughly matches
            # `huge_disk` (also used in the AS-VQSR VCF-gather job)
            worker_boot_disk_size=200,
            secondary_worker_boot_disk_size=200,
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
            gvcf_count=len(cohort.get_samples()),
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
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.load_vqsr import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                site_only_vcf_path=inputs.as_path(cohort, Vqsr, key='vcf'),
                vqsr_ht_path=self.expected_outputs(cohort),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[Combiner, SampleQC, Relatedness])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return get_workflow().prefix / 'frequencies.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.frequencies import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                vds_path=inputs.as_path(cohort, stage=Combiner),
                sample_qc_ht_path=inputs.as_path(cohort, stage=SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort, stage=Relatedness, key='relateds_to_drop'
                ),
                out_ht_path=self.expected_outputs(cohort),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
