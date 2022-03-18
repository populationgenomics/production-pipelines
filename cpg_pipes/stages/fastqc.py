"""
Stage that runs FastQC on alignment inputs.
"""

import logging

from cpg_pipes.storage import Path

from cpg_pipes.jobs import fastqc
from cpg_pipes.pipeline.dataset import Sample
from cpg_pipes.pipeline.pipeline import stage, SampleStage
from cpg_pipes.pipeline.exceptions import PipelineError
from cpg_pipes.pipeline.stage import StageInput, StageOutput

logger = logging.getLogger(__file__)


@stage
class FastQC(SampleStage):
    """
    Run FastQC on alignment inputs.
    """
    def expected_result(self, sample: Sample) -> dict[str, Path]:
        """
        Stage is expected to generate a FastQC HTML report, and a zip file for 
        parsing with MuiltiQC.
        """
        folder = sample.dataset.get_bucket() / 'qc'
        return {
            'html': folder / (sample.id + '_fastqc.html'),
            'zip': folder / (sample.id + '_fastqc.zip'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "fastqc" function implemented in the jobs module
        """
        if not sample.alignment_input:
            if self.skip_samples_with_missing_input:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(
                    f'No alignment input found for {sample.id}. '
                    f'Checked: Sequence entry and type=CRAM Analysis entry'
                )

        jobs = fastqc.fastqc(
            b=self.b,
            output_html_path=self.expected_result(sample)['html'],
            output_zip_path=self.expected_result(sample)['zip'],
            alignment_input=sample.alignment_input,
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=jobs
        )
