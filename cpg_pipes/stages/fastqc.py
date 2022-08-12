"""
Stage that runs FastQC on alignment inputs.
"""

import logging

from cpg_utils.config import get_config
from cpg_utils import Path

from cpg_pipes.jobs import fastqc
from cpg_pipes.targets import Sample
from cpg_pipes.pipeline import stage, SampleStage, StageInput, StageOutput
from cpg_pipes.filetypes import CramPath

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
        folder = sample.dataset.prefix() / 'qc'
        return {
            'html': folder / (sample.id + '_fastqc.html'),
            'zip': folder / (sample.id + '_fastqc.zip'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Only running FastQC if sequencing inputs are available.
        """
        sequencing_type = get_config()['workflow']['sequencing_type']
        alignment_input = sample.alignment_input_by_seq_type.get(sequencing_type)

        if isinstance(alignment_input, CramPath) and not alignment_input.is_bam:
            logger.info(
                f'FastQC input {sample} has CRAM inputs {alignment_input} '
                f'for type {sequencing_type}, skipping FASTQC'
            )
            return self.make_outputs(sample, skipped=True)

        if alignment_input is None or (
            get_config()['workflow'].get('check_inputs')
            and not alignment_input.exists()
        ):
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.error(f'No alignment inputs, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No FASTQ/BAM input found')

        jobs = fastqc.fastqc(
            b=self.b,
            output_html_path=self.expected_outputs(sample)['html'],
            output_zip_path=self.expected_outputs(sample)['zip'],
            alignment_input=alignment_input,
            job_attrs=self.get_job_attrs(sample),
            subsample=False,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
