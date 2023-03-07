"""
Stage to call SNVs in the mitochondrial genome of a single sample.

Reimplemented version of;
https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

"""
from typing import Any

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import mito, picard
from cpg_workflows.stages.align import Align
from cpg_workflows.targets import Sample
from cpg_workflows.utils import exists
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)


@stage(
    required_stages=Align,
    # analysis_type='vcf',
)
class AlignAndGenotypeMito(SampleStage):
    """
    Re-align and call SNVs in the mitochondrial genome of a single sample.

    This is a re-implementation of the broad pipeline as of 03/22 that was used for
    gnomad v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

    A second stage (not implemented here) will take the per sample outputs of this sage
    and merge them into a pseudo-joint-call set.

    Mitochondrial variant calling is a subtle art with potential for artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to avoid these
    issues, this stage has blind faith in the Broad pipeline and faithfully re-implements
    as much of the logic, tools and configuration as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sample cram.
        - Realignment of chrM reads to chrM reference using bwa.
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted"
            reference is the same sequence but with a different linearisation point. This
            is used to overcome mapping artifacts at the boundary (aka the "control
            region" of the mt genome).
        - Calling of variants from both crams using mutect2
        - Merging of the two call sets into a single vcf in normal chrM coordinate space.
        - Variant filtering part 1: Exclude black list of known problem sites and high
            alt allele counts.
        - haplocheckCLI is used to estimate contamination (using known mito haplotypes)
        - Variant filtering part 1: exclude variants with VAF below contamination
            estimate.
        - Export final vcf

    Mitochondrial reference indexes and filtering black lists were copied from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).

    Requires:
        Sample cram from Align stage.

    Outputs:
        out_vcf: the final filtered vcf for downstream use
        haplocheck_metrics: Metrics generated from the haplocheckCLI tool including an
            estimate of contamination and the predicted mitochondrial haplotype found.

    Configuration options:

    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        main = sample.dataset.prefix()
        analysis = sample.dataset.analysis_prefix()
        return {
            # Pipeline outputs
            'out_vcf': main / 'mito' / f'{sample.id}.mito.vcf',
            'haplocheck_metrics': analysis / 'mito' / f'{sample.id}.haplocheck.txt',
            'coverage': main / 'mito' / f'{sample.id}.coverage.tsv',
            # 'duplicate_metrics': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'coverage_metrics': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'theoretical_sensitivity_metrics': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'contamination_metrics': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'base_level_coverage_metrics': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'mean_coverage': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'median_coverage': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'major_haplogroup': sample.dataset.prefix() / f'{sample.id}.mito.foo',
            # 'contamination': sample.dataset.prefix() / f'{sample.id}.mito.foo',
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:

        # Mitochondrial specific bwa indexes.
        # QUESTION: Are these ok here or should they be defined centrally?
        mito_ref = get_batch().read_input_group(
            dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict',
            base='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta',
            amb='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb',
            ann='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann',
            bwt='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt',
            fai='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai',
            pac='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac',
            sa='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa',
        )
        shifted_mito_ref = get_batch().read_input_group(
            dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict',
            base='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta',
            amb='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb',
            ann='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann',
            bwt='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt',
            fai='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai',
            pac='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac',
            sa='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa',
            shift_back_chain='gs://cpg-common-main/references/hg38/v0/chrM/ShiftBack.chain',
        )
        intervals = get_batch().read_input_group(
            control_region_shifted='gs://cpg-common-main/references/hg38/v0/chrM/control_region_shifted.chrM.interval_list',
            non_control_region='gs://cpg-common-main/references/hg38/v0/chrM/non_control_region.chrM.interval_list',
        )


        jobs = []

        # Extract reads mapped to chrM
        cram_path = inputs.as_path(sample, Align, 'cram')
        crai_path = inputs.as_path(sample, Align, 'crai')
        subset_bam_j = mito.subset_cram_to_chrM(
            b=get_batch(),
            cram_path=CramPath(cram_path, crai_path),
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(subset_bam_j)

        # Align extracted reads to chrM using bwa
        realign_j = mito.mito_realign(
            b=get_batch(),
            sample_id=sample.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=mito_ref,
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(realign_j)

        # Mark duplicates (non-shifted)
        realign_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=realign_j.output_cram,
            fasta_reference=mito_ref,
        )
        jobs.append(realign_mkdup_j)

        # Align extracted reads to "shifted" chrM using bwa
        shifted_realign_j = mito.mito_realign(
            b=get_batch(),
            sample_id=sample.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(shifted_realign_j)

        # Mark duplicates (shifted)
        shifted_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=shifted_realign_j.output_cram,
            fasta_reference=mito_ref,
        )
        jobs.append(shifted_mkdup_j)

        # Call variants on WT genome
        call_j = mito.mito_mutect2(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            reference=mito_ref,
            region='chrM:576-16024',  # Exclude the control region.
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(call_j)

        # Call variants in ONLY the control region using the shifted reference
        shifted_call_j = mito.mito_mutect2(
            b=get_batch(),
            cram=shifted_mkdup_j.output_cram,
            reference=shifted_mito_ref,
            region='chrM:8025-9144',  # Only call inside the control region.
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(shifted_call_j)

        # Merge the wt and shifted VCFs
        merge_j = mito.liftover_and_combine_vcfs(
            b=get_batch(),
            vcf=call_j.output_vcf,
            shifted_vcf=shifted_call_j.output_vcf,
            reference=mito_ref,
            shift_back_chain=shifted_mito_ref.shift_back_chain,
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(merge_j)

        # Merge the mutect stats output files (needed for filtering)
        merge_stats_J = mito.merge_mutect_stats(
            b=get_batch(),
            first_stats_file=call_j.output_vcf['vcf.gz.stats'],
            second_stats_file=shifted_call_j.output_vcf['vcf.gz.stats'],
        )
        jobs.append(merge_stats_J)

        # Initial round of filtering to exclude blacklist and high alt alleles
        initial_filter_j = mito.filter_variants(
            b=get_batch(),
            vcf=merge_j.output_vcf,
            reference=mito_ref,
            merged_mutect_stats=merge_stats_J.combined_stats,
            # alt_allele and vaf config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=4,
            min_allele_fraction=0,
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(initial_filter_j)

        # SplitMultiAllelics AND remove non-passing sites
        # Output is only used for input to haplocheck
        split_multiallelics_j = mito.split_multi_allelics(
            b=get_batch(),
            vcf=initial_filter_j.output_vcf,
            reference=mito_ref,
            remove_non_pass_sites=True
        )
        jobs.append(split_multiallelics_j)

        # Use mito reads to identify level of contamination
        get_contamination_j = mito.get_contamination(
            b=get_batch(),
            vcf=split_multiallelics_j.output_vcf,
            haplocheck_output=self.expected_outputs(sample)['haplocheck_metrics'],
        )
        jobs.append(get_contamination_j)

        # Filter round 2 - remove variants with VAF below estimated contamination
        second_filter_j = mito.filter_variants(
            b=get_batch(),
            vcf=initial_filter_j.output_vcf,
            reference=mito_ref,
            merged_mutect_stats=merge_stats_J.combined_stats,
            # alt_allele and vaf config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=4,
            min_allele_fraction=0,
            contamination_estimate=get_contamination_j.max_contamination,
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(second_filter_j)

        # CoverageAtEveryBase
        # This is a temporary task to handle "joint calling" until Mutect2 can produce a
        # GVCF. This proivdes coverage at each base so low coverage sites can be
        # considered ./. rather than 0/0.
        # Generate final output vcf
        non_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=realign_j.output_cram,
            intervals_list=intervals.non_control_region,
            reference=mito_ref,
        )
        jobs.append(non_control_region_coverage_j)

        shifted_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=realign_j.output_cram,
            intervals_list=intervals.control_region_shifted,
            reference=mito_ref,
        )
        jobs.append(shifted_control_region_coverage_j)

        # Merge coverage stats
        merge_coverage_j = mito.merge_coverage(
            b=get_batch(),
            non_cr_coverage=non_control_region_coverage_j.per_base_coverage,
            shifted_cr_coverage=shifted_control_region_coverage_j.per_base_coverage,
            merged_coverage=self.expected_outputs(sample)['coverage']
        )
        jobs.append(merge_coverage_j)

        # Generate final output vcf
        split_multiallelics_j = mito.split_multi_allelics(
            b=get_batch(),
            vcf=second_filter_j.output_vcf,
            reference=mito_ref,
            remove_non_pass_sites=False
        )
        jobs.append(split_multiallelics_j)

        # Write the final vcf to the bucket
        get_batch().write_output(j.split_multiallelics_j.output_vcf, str(out_vcf))

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
