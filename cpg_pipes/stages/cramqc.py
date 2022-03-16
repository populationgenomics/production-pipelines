"""
Stages that perform alignment QC on CRAM files.
"""

import logging

from cpg_pipes.jobs.cram_qc import samtools_stats, verify_bamid, picard_wgs_metrics
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput

logger = logging.getLogger(__file__)


@stage
class SamtoolsStats(SampleStage):
    """
    Alignment QC using samtools stats.
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate one QC file to be parsed with MultiQC:
        * Samtools stats file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L652-L654
        """
        return sample.dataset.get_bucket() / 'qc' / (sample.id + '_samtools_stats.txt')

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = sample.get_cram_path()

        j = samtools_stats(
            b=self.pipe.b,
            cram_path=cram_path,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(sample, data=self.expected_result(sample), jobs=[j])


@stage
class PicardWgsMetrics(SampleStage):
    """
    Alignment QC using picard tools wgs_metrics.
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate one QC file to be parsed with MultiQC.
        * Picard file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L539-L541
        """
        return (
            sample.dataset.get_bucket() / 'qc' / (sample.id + '_picard_wgs_metrics.csv')
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = sample.get_cram_path()

        j = picard_wgs_metrics(
            b=self.pipe.b,
            cram_path=cram_path,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(sample, data=self.expected_result(sample), jobs=[j])


@stage
class VerifyBamId(SampleStage):
    """
    Check for contamination using VerifyBAMID SelfSM.
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate one QC file to be parsed with MultiQC.
        * VerifyBAMID file has to have *.selfSM ending:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L783-L784
        """
        return (
            sample.dataset.get_bucket() / 'qc' / (sample.id + '_verify_bamid.selfSM')
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = sample.get_cram_path()

        j = verify_bamid(
            b=self.pipe.b,
            cram_path=cram_path,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(sample, data=self.expected_result(sample), jobs=[j])
