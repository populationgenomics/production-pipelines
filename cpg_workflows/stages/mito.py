"""
Stage to call SNVs in the mitochondrial genome of a single sequencing_group.

Reimplemented version of;
https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

"""
from typing import Any

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import mito, picard
from cpg_workflows.stages.align import Align
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroupStage,
    SequencingGroup,
)

MITO_REF = {
    'dict': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict',
    'base': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta',
    'amb': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb',
    'ann': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann',
    'bwt': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt',
    'fai': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai',
    'pac': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac',
    'sa': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa',
}

SHIFTED_MITO_REF = {
    'dict': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict',
    'base': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta',
    'amb': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb',
    'ann': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann',
    'bwt': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt',
    'fai': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai',
    'pac': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac',
    'sa': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa',
    'shift_back_chain': 'gs://cpg-common-main/references/hg38/v0/chrM/ShiftBack.chain',
}

CONTROL_REGION_INTERVALS = {
    'control_region_shifted': 'gs://cpg-common-main/references/hg38/v0/chrM/control_region_shifted.chrM.interval_list',
    'non_control_region': 'gs://cpg-common-main/references/hg38/v0/chrM/non_control_region.chrM.interval_list',
}


@stage(
    required_stages=Align,
    # TODO: add suitable analysis types
)
class RealignMito(SequencingGroupStage):
    """
    Re-align mitochondrial genome of a single sequencing_group.

    This is a re-implementation of the broad pipeline (as of 03/22) used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

    Mitochondrial variant calling can be subject to subtle artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to minimize
    these issues, this pipeline attempts to re-implement as much of the logic, tools and
    configuration as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sequencing_group cram.
        - Realignment of chrM reads to chrM reference using bwa.
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted"
            reference is the same sequence but with a different linearisation point. This
            is used to overcome mapping artifacts at the boundary (aka the "control
            region" of the mt genome).
        - Generation of overall and  base level alignment statistics.

    Mitochondrial reference indexes and lift over chains were sourced from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).

    Requires:
        sequencing_group cram from Align stage.

    Outputs:
        non_shifted_cram: Sorted cram with duplicates marked aligned to the standard
            hg38 chrM.
        shifted_cram: Sorted cram with duplicates marked aligned to a copy of chrM
            with the linearisation location "shifted"
        base_level_coverage_metrics: per base coverage needed to differentiate between
            0/0 and ./. genotypes until mutect can generate a gVCF.
        coverage_metrics: CollectWgsMetrics output.
        theoretical_sensitivity_metrics: CollectWgsMetrics output.

    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        analysis = sequencing_group.dataset.analysis_prefix()
        return {
            'non_shifted_cram': main / 'mito' / f'{sequencing_group.id}.mito.cram',
            'shifted_cram': main / 'mito' / f'{sequencing_group.id}.shifted_mito.cram',
            'base_level_coverage_metrics': main
            / 'mito'
            / f'{sequencing_group.id}.base_level_coverage.tsv',
            'coverage_metrics': analysis
            / 'mito'
            / f'{sequencing_group.id}.coverage_metrics.txt',
            'coverage_mean': analysis
            / 'mito'
            / f'{sequencing_group.id}.coverage_mean.txt',
            'coverage_median': analysis
            / 'mito'
            / f'{sequencing_group.id}.coverage_median.txt',
            'theoretical_sensitivity_metrics': analysis
            / 'mito'
            / f'{sequencing_group.id}.theoretical_sensitivity.txt',
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # Mitochondrial specific reference files.
        mito_ref = get_batch().read_input_group(**MITO_REF)
        shifted_mito_ref = get_batch().read_input_group(**SHIFTED_MITO_REF)
        intervals = get_batch().read_input_group(**CONTROL_REGION_INTERVALS)

        jobs = []

        # Extract reads mapped to chrM
        cram_path = inputs.as_path(sequencing_group, Align, 'cram')
        crai_path = inputs.as_path(sequencing_group, Align, 'crai')
        subset_bam_j = mito.subset_cram_to_chrM(
            b=get_batch(),
            cram_path=CramPath(cram_path, crai_path),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(subset_bam_j)
        assert isinstance(subset_bam_j.output_bam, hb.ResourceGroup)

        # Align extracted reads to chrM using bwa
        realign_j = mito.mito_realign(
            b=get_batch(),
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(realign_j)
        assert isinstance(realign_j.output_cram, hb.ResourceFile)

        # Mark duplicates (non-shifted)
        realign_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=realign_j.output_cram,
            fasta_reference=mito_ref,
            output_path=self.expected_outputs(sequencing_group)['non_shifted_cram'],
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=True,
        )
        jobs.append(realign_mkdup_j)
        assert isinstance(realign_mkdup_j, Job)
        assert isinstance(realign_mkdup_j.output_cram, hb.ResourceGroup)

        # Align extracted reads to "shifted" chrM using bwa
        shifted_realign_j = mito.mito_realign(
            b=get_batch(),
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_realign_j)
        assert isinstance(shifted_realign_j.output_cram, hb.ResourceFile)

        # Mark duplicates (shifted)
        shifted_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=shifted_realign_j.output_cram,
            fasta_reference=shifted_mito_ref,
            output_path=self.expected_outputs(sequencing_group)['shifted_cram'],
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=True,
        )
        jobs.append(shifted_mkdup_j)
        assert isinstance(shifted_mkdup_j, Job)
        assert isinstance(shifted_mkdup_j.output_cram, hb.ResourceGroup)

        # Collect coverage metrics (only on non-shifted)
        coverage_metrics_J = mito.collect_coverage_metrics(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            reference=mito_ref,
            metrics=self.expected_outputs(sequencing_group)['coverage_metrics'],
            theoretical_sensitivity=self.expected_outputs(sequencing_group)[
                'theoretical_sensitivity_metrics'
            ],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(coverage_metrics_J)
        assert isinstance(coverage_metrics_J.metrics, hb.ResourceFile)

        # Extract mean and median coverage (only on non-shifted)
        extract_coverage_mean_j = mito.extract_coverage_mean(
            b=get_batch(),
            metrics=coverage_metrics_J.metrics,
            mean_path=self.expected_outputs(sequencing_group)['coverage_mean'],
            median_path=self.expected_outputs(sequencing_group)['coverage_median'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(extract_coverage_mean_j)

        # CoverageAtEveryBase
        # Provides coverage at each base so low coverage sites can be
        # considered ./. rather than 0/0.
        non_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            intervals_list=intervals.non_control_region,
            reference=mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(non_control_region_coverage_j)
        assert isinstance(
            non_control_region_coverage_j.per_base_coverage, hb.ResourceFile
        )

        shifted_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=shifted_mkdup_j.output_cram,
            intervals_list=intervals.control_region_shifted,
            reference=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_control_region_coverage_j)
        assert isinstance(
            shifted_control_region_coverage_j.per_base_coverage, hb.ResourceFile
        )

        # Merge coverage stats
        merge_coverage_j = mito.merge_coverage(
            b=get_batch(),
            non_cr_coverage=non_control_region_coverage_j.per_base_coverage,
            shifted_cr_coverage=shifted_control_region_coverage_j.per_base_coverage,
            merged_coverage=self.expected_outputs(sequencing_group)[
                'base_level_coverage_metrics'
            ],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_coverage_j)

        return self.make_outputs(
            sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs
        )
