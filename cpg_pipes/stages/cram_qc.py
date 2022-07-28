"""
Stages that perform alignment QC on CRAM files.
"""

import logging

from cpg_utils import Path

from ..jobs.cram_qc import samtools_stats, verify_bamid, picard_wgs_metrics
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from ..targets import Sample

logger = logging.getLogger(__file__)


@stage
class SamtoolsStats(SampleStage):
    """
    Alignment QC using samtools stats.
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Expected to generate one QC file to be parsed with MultiQC:
        * Samtools stats file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L652-L654
        """
        return sample.dataset.prefix() / 'qc' / (sample.id + '_samtools_stats.txt')

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Call a function from the `jobs` module.
        """
        cram_path = sample.get_cram_path()

        j = samtools_stats(
            b=self.b,
            cram_path=cram_path,
            output_path=self.expected_outputs(sample),
            job_attrs=self.get_job_attrs(sample),
            sequencing_type=self.cohort.sequencing_type,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[j])


@stage
class PicardWgsMetrics(SampleStage):
    """
    Alignment QC using picard tools wgs_metrics.
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Expected to generate one QC file to be parsed with MultiQC.
        * Picard file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L539-L541
        """
        return sample.dataset.prefix() / 'qc' / (sample.id + '_picard_wgs_metrics.csv')

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Call a function from the `jobs` module.
        """
        cram_path = sample.get_cram_path()

        j = picard_wgs_metrics(
            b=self.b,
            cram_path=cram_path,
            output_path=self.expected_outputs(sample),
            job_attrs=self.get_job_attrs(sample),
            sequencing_type=self.cohort.sequencing_type,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[j])


@stage
class VerifyBamId(SampleStage):
    """
    Check for contamination using VerifyBAMID SelfSM.
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Expected to generate one QC file to be parsed with MultiQC.
        * VerifyBAMID file has to have *.selfSM ending:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L783-L784
        """
        return sample.dataset.prefix() / 'qc' / (sample.id + '_verify_bamid.selfSM')

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Call a function from the `jobs` module.
        """
        cram_path = sample.get_cram_path()

        j = verify_bamid(
            b=self.b,
            cram_path=cram_path,
            output_path=self.expected_outputs(sample),
            job_attrs=self.get_job_attrs(sample),
            sequencing_type=self.cohort.sequencing_type,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[j])
