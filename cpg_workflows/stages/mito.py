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

    This is a re-implementation of the broad pipeline as of 03/22 that was used for gnomad v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

    A second stage (not implemented here) will take the per sample outputs of this sage and merge them
    into a psudo-joint-call set.

    Mitochondrial variant calling is a subtle art with potential for many artifacts resulting from mapping
    errors and other complexities. To avoid introducing these issues this stage faithfully re-implements
    as much of the logic, tools and configuration used in the Broad pipeline as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sample cram.
        - Realignment of chrM reads to chrM reference using bwa
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted" reference is the
            same sequence but with a different linearisation point. This is used to overcome mapping artifacts
            at the boundary.
        - Calling of variants from both crams using mutect2
        - Merging of the two call sets into a single vcf in normal chrM coordinate space.
        - Variant filtering part 1: Exclude black list of known problem sites and high alt allele counts.
        - Run haplocheck to estimate contamination (using known mito haplotypes)

    Expected outputs include:
     - out_vcf: the final filtered vcf for downstream use
     - haplocheck_metrics: Metrics generated from the haplocheckCLI tool including an estimate of contamination and the predicted
        mitochondrial haplotype found.
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        # Listing all outputs from Broad WDL here commenting what we chose not to generate
        return {
            'mito_subset_bam': sample.dataset.prefix() / 'mito' / f'{sample.id}.mito_subset.bam',
            'mito_realigned_cram': sample.dataset.prefix() / 'mito' / f'{sample.id}.mito_realign.cram',
            'mito_shifted_cram': sample.dataset.prefix() / 'mito' / f'{sample.id}.mito_shifted.cram',
            'out_vcf': sample.dataset.prefix() / 'mito' / f'{sample.id}.mito.vcf',
            'haplocheck_metrics': sample.dataset.prefix() / 'mito' / f'{sample.id}.haplocheck.txt',
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
        cram_path = inputs.as_path(sample, Align, 'cram')
        crai_path = inputs.as_path(sample, Align, 'crai')

        jobs = []

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

        # Extract reads mapped to chrM
        j, subset_bam = mito.subset_cram_to_chrM(
            b=get_batch(),
            cram_path=CramPath(cram_path, crai_path),
            output_bam_path=self.expected_outputs(sample)['mito_subset_bam'],
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(j)

        # Align to chrM genome using bwa
        jobs += mito.mito_realign(
            b=get_batch(),
            sample_id=sample.id,
            input_bam=subset_bam,
            output_cram_path=self.expected_outputs(sample)['mito_realigned_cram'],
            mito_ref=mito_ref,
            job_attrs=self.get_job_attrs(sample),
        )

        # Align to chrM genome with offset linearisation site to account for
        # alignment boundry effects at start and end of reff genome
        jobs += mito.mito_realign(
            b=get_batch(),
            sample_id=sample.id,
            input_bam=subset_bam,
            output_cram_path=self.expected_outputs(sample)['mito_shifted_cram'],
            mito_ref=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sample),
        )

        # Call variants
        genotype_jobs = mito.genotype_mito(
            b=get_batch(),
            cram_path=self.expected_outputs(sample)['mito_realigned_cram'],
            shifted_cram_path=self.expected_outputs(sample)['mito_shifted_cram'],
            output_vcf_path=self.expected_outputs(sample)['out_vcf'],
            mito_reff=mito_ref,
            shifted_mito_reff=shifted_mito_ref,
            haplocheck_output=self.expected_outputs(sample)['haplocheck_output']
        )
        for job in genotype_jobs:
            if job:
                job.depends_on(*[j for j in jobs if j])
        jobs += genotype_jobs

        # InitialFilter
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
