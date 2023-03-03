"""
Stage to call mito genome variants

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
    Re-align and call mitochondrial genome using bwa and mutect2
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        # Listing all outputs from Broad WDL here commenting what we chose not to generate
        return {
    'mito_subset_bam': sample.dataset.prefix() / f'{sample.id}.mito_subset.bam',
    'mito_realigned_cram': sample.dataset.prefix() / f'{sample.id}.mito_realign.cram',
    'mito_realigned_crai': sample.dataset.prefix() / f'{sample.id}.mito_realign.cram.crai',
    'mito_shifted_cram': sample.dataset.prefix() / f'{sample.id}.mito_shifted.cram',
    'mito_shifted_crai': sample.dataset.prefix() / f'{sample.id}.mito_shifted.cram.crai',
    'out_vcf': sample.dataset.prefix() /  f'{sample.id}.mito.vcf',
    # 'out_vcf_index': sample.dataset.prefix() / f'{sample.id}.mito.foo',
    # 'split_vcf': sample.dataset.prefix() / f'{sample.id}.mito.foo',
    # 'split_vcf_index': sample.dataset.prefix() / f'{sample.id}.mito.foo',
    # 'input_vcf_for_haplochecker': sample.dataset.prefix() / f'{sample.id}.mito.foo',
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
            input_bam=subset_bam,
            output_cram_path=self.expected_outputs(sample)['mito_realigned_cram'],
            shifted=False,
            job_attrs=self.get_job_attrs(sample),
        )

        # Align to chrM genome with offset linearisation site to account for
        # alignment boundry effects at start and end of reff genome
        jobs += mito.mito_realign(
            b=get_batch(),
            input_bam=subset_bam,
            output_cram_path=self.expected_outputs(sample)['mito_shifted_cram'],
            shifted=True,
            job_attrs=self.get_job_attrs(sample),
        )

        # Call variants




        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
