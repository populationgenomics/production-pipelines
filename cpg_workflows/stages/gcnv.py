"""
Stages that implement GATK-gCNV.
"""
from typing import Any

from cpg_utils import Path
from cpg_utils.config import get_config, try_get_ar_guid, AR_GUID_NAME
from cpg_workflows.jobs import gcnv
from cpg_workflows.targets import SequencingGroup, Cohort
from cpg_workflows.workflow import (
    stage,
    CohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
)

from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    make_combined_ped,
    get_images,
    get_references,
    queue_annotate_sv_jobs,
    _gcnv_annotated_meta,
)

from cpg_workflows.batch import get_batch


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
            'intervals': self.prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'intervals_index': self.prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'segments': self.prefix / f'{seqgroup.id}.segments.vcf.gz',
            'segments_index': self.prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'ratios': self.prefix / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(
        self, seqgroup: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
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


@stage(required_stages=[GermlineCNVCalls])
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
        gcnv_vcfs = inputs.as_dict_by_target(GermlineCNVCalls)
        all_vcfs = [
            str(gcnv_vcfs[sgid]['intervals'])
            for sgid in cohort.get_sequencing_group_ids()
        ]

        job_or_none = gcnv.merge_calls(
            get_batch(),
            sg_vcfs=all_vcfs,
            job_attrs=self.get_job_attrs(cohort),
            output_path=outputs['combined_calls'],
        )
        return self.make_outputs(cohort, data=outputs, jobs=job_or_none)


@stage(required_stages=FastCombineGCNVs)
class TranslategCNVToGATKVCF(CohortStage):
    """
    attempt at fixing the gCNV workflow - pass the VCF through
    the to-GATK format VCF parser. This operates the Cromwell
    workflow from GATK-SV, and may not work at all...
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Args:
            cohort (Cohort): cohort of all samples (across several sub-cohort batches)
            inputs (StageInput): access to prior inputs
        """

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'vcf': inputs.as_dict(cohort, FastCombineGCNVs)['combined_calls'],
            'ped_file': make_combined_ped(cohort, self.prefix),
        }
        input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(cohort)

        billing_labels = {
            'stage': 'formatvcfforgatk',
            AR_GUID_NAME: try_get_ar_guid(),
        }

        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name='FormatVcfForGatk',
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(
    required_stages=TranslategCNVToGATKVCF,
    analysis_type='sv',
    analysis_keys=['annotated_vcf'],
    update_analysis_meta=_gcnv_annotated_meta,
)
class AnnotateVcfSV(CohortStage):
    """
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
            input_vcf=inputs.as_dict(cohort, TranslategCNVToGATKVCF)[
                'gatk_formatted_vcf'
            ],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_out, jobs=job_or_none)
