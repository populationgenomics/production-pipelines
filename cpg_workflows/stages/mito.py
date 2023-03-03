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
    'mito_subset_cram': sample.dataset.prefix() / f'{sample.id}.mito_subset.cram',
    'mito_subset_crai': sample.dataset.prefix() / f'{sample.id}.mito_subset.cram.crai',
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
        # Generate minibam
        # cpg_workflows/jobs/align.py: extract_fastq
        # sort_cmd
        # finalise_alignment (mark dups)

        # map twice
        # collect metrics

        j = mito.subset_cram_to_chrM(
            b=get_batch(),
            cram_path=CramPath(cram_path, crai_path),
            mito_subset_cram=self.expected_outputs(sample)['mito_subset_cram'],
            mito_subset_crai=self.expected_outputs(sample)['mito_subset_crai'],
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(j)

        jobs += mito.mito_realign(
            b=get_batch(),
            cram_path=CramPath(self.expected_outputs(sample)['mito_subset_cram'], self.expected_outputs(sample)['mito_subset_crai']),
            mito_aligned_cram=self.expected_outputs(sample)['mito_realigned_cram'],
            mito_aligned_crai=self.expected_outputs(sample)['mito_realigned_crai'],
            shifted=False,
            job_attrs=self.get_job_attrs(sample),
        )

        jobs += mito.mito_realign(
            b=get_batch(),
            cram_path=CramPath(self.expected_outputs(sample)['mito_subset_cram'], self.expected_outputs(sample)['mito_subset_crai']),
            mito_aligned_cram=self.expected_outputs(sample)['mito_shifted_cram'],
            mito_aligned_crai=self.expected_outputs(sample)['mito_shifted_crai'],
            shifted=True,
            job_attrs=self.get_job_attrs(sample),
        )

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
