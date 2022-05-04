"""
Stage that runs FastQC on alignment inputs.
"""

import logging

from .. import Path, types
from ..jobs import fastqc
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
        folder = sample.dataset.get_bucket() / 'qc'
        return {
            'html': folder / (sample.id + '_fastqc.html'),
            'zip': folder / (sample.id + '_fastqc.zip'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Only running FastQC if sequencing inputs are available.
        """
        alignment_input = sample.alignment_input
        if not alignment_input or (
            self.check_inputs and not types.alignment_input_exists(alignment_input)
        ):
            if (
                self.pipeline_config.get('fastqc_on_realigned_cram', False)
                and sample.get_cram_path().exists()
            ):
                logger.info(
                    f'No alignment inputs found for sample {sample.id}, using '
                    f'aligned CRAM instead'
                )
                alignment_input = sample.get_cram_path()
            else:
                logger.warning(f'No alignment inputs, skipping sample {sample.id}')
                return None

        jobs = fastqc.fastqc(
            b=self.b,
            output_html_path=self.expected_outputs(sample)['html'],
            output_zip_path=self.expected_outputs(sample)['zip'],
            alignment_input=alignment_input,
            refs=self.refs,
            job_attrs=self.get_job_attrs(sample),
            subsample=False,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
