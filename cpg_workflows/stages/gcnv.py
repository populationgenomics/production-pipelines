"""
Stages that implement GATK-gCNV.
"""

from cpg_utils import Path
from cpg_utils.config import get_config, try_get_ar_guid, AR_GUID_NAME
from cpg_utils.hail_batch import get_batch, image_path
from cpg_workflows.jobs import gcnv
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    get_images,
    get_references,
    queue_annotate_sv_jobs,
)
from cpg_workflows.targets import SequencingGroup, Cohort
from cpg_workflows.workflow import (
    stage,
    CohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
)


def _gcnv_annotated_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Callable, adds custom analysis object meta attribute
    """
    return {'type': 'gCNV-annotated'}


@stage
class PrepareIntervals(CohortStage):
    """
    Interval preparation steps that don't require the sample read counts:
    PreprocessIntervals and AnnotateIntervals.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'preprocessed': self.prefix / f'{cohort.name}.preprocessed.interval_list',
            'annotated': self.prefix / f'{cohort.name}.annotated.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.prepare_intervals(
            get_batch(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=PrepareIntervals)
class CollectReadCounts(SequencingGroupStage):
    """
    Per-sample stage that runs CollectReadCounts to produce .counts.tsv.gz files.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'counts': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.counts.tsv.gz',
            'index': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.counts.tsv.gz.tbi',
        }

    def queue_jobs(
        self, seqgroup: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        if seqgroup.cram is None:
            raise ValueError(f'No CRAM file found for {seqgroup}')

        jobs = gcnv.collect_read_counts(
            get_batch(),
            inputs.as_path(seqgroup.dataset, PrepareIntervals, 'preprocessed'),
            seqgroup.cram,
            self.get_job_attrs(seqgroup),
            seqgroup.dataset.prefix() / 'gcnv' / seqgroup.id,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=[PrepareIntervals, CollectReadCounts])
class DeterminePloidy(CohortStage):
    """
    The non-sharded cohort-wide gCNV steps after read counts have been collected:
    FilterIntervals and DetermineGermlineContigPloidy. These outputs represent
    intermediate results for the cohort as a whole, so are written to tmp_prefix.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'filtered': self.tmp_prefix / f'{cohort.name}.filtered.interval_list',
            'calls': self.tmp_prefix / f'{cohort.name}-ploidy-calls.tar.gz',
            'model': self.tmp_prefix / f'{cohort.name}-ploidy-model.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.filter_and_determine_ploidy(
            get_batch(),
            get_config()['workflow'].get('ploidy_priors'),
            inputs.as_path(cohort, PrepareIntervals, 'preprocessed'),
            inputs.as_path(cohort, PrepareIntervals, 'annotated'),
            inputs.as_path_by_target(CollectReadCounts, 'counts').values(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[PrepareIntervals, CollectReadCounts, DeterminePloidy])
class GermlineCNV(CohortStage):
    """
    The cohort-wide GermlineCNVCaller step, sharded across genome regions.
    This is separate from the DeterminePloidy stage so that the GermlineCNVCalls
    stage can pick out this stage's sharded inputs easily.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            name: self.tmp_prefix / f'{name}.tar.gz' for name in gcnv.shard_basenames()
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.shard_gcnv(
            get_batch(),
            inputs.as_path(cohort, PrepareIntervals, 'annotated'),
            inputs.as_path(cohort, DeterminePloidy, 'filtered'),
            inputs.as_path(cohort, DeterminePloidy, 'calls'),
            inputs.as_path_by_target(CollectReadCounts, 'counts').values(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[DeterminePloidy, GermlineCNV])
class GermlineCNVCalls(SequencingGroupStage):
    """
    Produces final individual VCF results by running PostprocessGermlineCNVCalls.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'intervals': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.intervals.vcf.gz',
            'intervals_index': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'segments': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.segments.vcf.gz',
            'segments_index': seqgroup.dataset.prefix()
            / 'gcnv'
            / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'ratios': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(seqgroup)

        jobs = gcnv.postprocess_calls(
            get_batch(),
            inputs.as_path(seqgroup.dataset, DeterminePloidy, 'calls'),
            inputs.as_dict(seqgroup.dataset, GermlineCNV),
            # FIXME get the sample index via sample_name.txt files instead
            seqgroup.dataset.get_sequencing_group_ids().index(seqgroup.id),
            self.get_job_attrs(seqgroup),
            output_prefix=str(self.prefix / seqgroup.id),
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=GermlineCNVCalls)
class PrepareVcfsForMerge(SequencingGroupStage):
    """
    Fixes the GT header in the VCFs produced by gCNV
    During experimentation these were found to be Integer instead of String
    This will be fixed in a future release of GATK
    Also splits multiallelic loci (all sites are called as Dup/Del GTs)
    Each site needs to have only one Alt allele for the annotation step
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'fixed_intervals': self.prefix / f'{seqgroup.id}.fixed_intervals.vcf.bgz',
            'fixed_intervals_index': self.prefix
            / f'{seqgroup.id}.fixed_intervals.vcf.bgz.tbi',
        }

    def queue_jobs(
        self, seqgroup: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        jobs = gcnv.fix_intervals_vcf(
            get_batch(),
            inputs.as_path(seqgroup, GermlineCNVCalls, 'intervals'),
            self.get_job_attrs(seqgroup),
            output_path=outputs['fixed_intervals'],
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=PrepareVcfsForMerge)
class FastCombineGCNVs(CohortStage):
    """
    Produces final multi-sample VCF results by running a merge
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'combined_calls': self.prefix / 'gcnv_joint_call.vcf.bgz',
            'combined_calls_index': self.prefix / 'gcnv_joint_call.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        # do a slapdash bcftools merge on all input files...
        gcnv_vcfs = inputs.as_dict_by_target(PrepareVcfsForMerge)
        all_vcfs = [
            str(gcnv_vcfs[sgid]['fixed_intervals'])
            for sgid in cohort.get_sequencing_group_ids()
        ]

        pipeline_image = get_images(['sv_pipeline_docker'])['sv_pipeline_docker']

        job_or_none = gcnv.merge_calls(
            get_batch(),
            sg_vcfs=all_vcfs,
            docker_image=pipeline_image,
            job_attrs=self.get_job_attrs(cohort),
            output_path=outputs['combined_calls'],
        )
        return self.make_outputs(cohort, data=outputs, jobs=job_or_none)


@stage(
    required_stages=FastCombineGCNVs,
    analysis_type='sv',
    analysis_keys=['annotated_vcf'],
    update_analysis_meta=_gcnv_annotated_meta,
)
class AnnotateCNV(CohortStage):
    """
    Smaller, direct annotation using SvAnnotate
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    This is a full clone of the GATK-SV pipeline Cromwell stage, but use on a slightly
    different output. Trying to work out the best way to handle this through inheritance

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        return {
            'annotated_vcf': self.prefix / 'unfiltered_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'unfiltered_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        expected_out = self.expected_outputs(cohort)

        billing_labels = {
            'stage': self.name.lower(),
            AR_GUID_NAME: try_get_ar_guid(),
        }

        job_or_none = queue_annotate_sv_jobs(
            batch=get_batch(),
            cohort=cohort,
            cohort_prefix=self.prefix,
            input_vcf=inputs.as_dict(cohort, FastCombineGCNVs)['combined_calls'],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_out, jobs=job_or_none)


@stage(required_stages=AnnotateCNV)
class AnnotateCNVVcfWithStrvctvre(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.prefix / 'cnv_strvctvre_annotated.vcf.bgz',
            'strvctvre_vcf_index': self.prefix / 'cnv_strvctvre_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        strv_job = get_batch().new_job(
            'StrVCTVRE', self.get_job_attrs() | {'tool': 'strvctvre'}
        )

        strv_job.image(image_path('strvctvre'))
        strv_job.storage('20Gi')

        strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']
        phylop_in_batch = get_batch().read_input(strvctvre_phylop)

        input_dict = inputs.as_dict(cohort, AnnotateCNV)
        expected_d = self.expected_outputs(cohort)

        # read vcf and index into the batch
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']

        strv_job.declare_resource_group(
            output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )

        # run strvctvre
        strv_job.command(
            f'python StrVCTVRE.py '
            f'-i {input_vcf} '
            f'-o {strv_job.output_vcf["vcf.gz"]} '
            f'-f vcf '
            f'-p {phylop_in_batch}'
        )
        strv_job.command(f'tabix {strv_job.output_vcf["vcf.gz"]}')

        get_batch().write_output(
            strv_job.output_vcf, str(expected_d['strvctvre_vcf']).replace('.vcf.gz', '')
        )
        return self.make_outputs(cohort, data=expected_d, jobs=strv_job)
