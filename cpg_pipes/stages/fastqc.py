"""
Stage that runs FastQC on alignment inputs.
"""

import logging

from .. import Path, types
from ..jobs import fastqc
from ..jobs.align import process_alignment_input
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput

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
        folder = sample.dataset.path() / 'qc'
        return {
            'html': folder / (sample.id + '_fastqc.html'),
            'zip': folder / (sample.id + '_fastqc.zip'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Only running FastQC if sequencing inputs are available.
        """
        seq_type, alignment_input = process_alignment_input(
            sample, 
            seq_type=self.pipeline_config.get('sequencing_type'),
            realign_cram_ver=self.pipeline_config.get('realign_from_cram_version'),
        )

        if alignment_input is None or (
            self.check_inputs and not alignment_input.exists()
        ):
            if self.skip_samples_with_missing_input:
                logger.error(f'No alignment inputs, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, error_msg=f'No alignment input found for {sample.id}'
                )

        jobs = fastqc.fastqc(
            b=self.b,
            output_html_path=self.expected_outputs(sample)['html'],
            output_zip_path=self.expected_outputs(sample)['zip'],
            alignment_input=alignment_input,
            refs=self.refs,
            images=self.images,
            job_attrs=self.get_job_attrs(sample),
            subsample=False,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
