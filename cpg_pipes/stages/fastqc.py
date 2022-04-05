"""
Stage that runs FastQC on alignment inputs.
"""

import logging

from .. import Path
from ..jobs import fastqc
from cpg_pipes.targets import Sample
from ..pipeline import stage, SampleStage, PipelineError, StageInput, StageOutput

logger = logging.getLogger(__file__)


@stage
class FastQC(SampleStage):
    """
    Run FastQC on alignment inputs.
    """
    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
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
            output_html_path=self.expected_outputs(sample)['html'],
            output_zip_path=self.expected_outputs(sample)['zip'],
            alignment_input=sample.alignment_input,
            refs=self.refs,
            job_attrs=sample.get_job_attrs(),
        )
        return self.make_outputs(
            sample, 
            data=self.expected_outputs(sample), 
            jobs=jobs
        )
