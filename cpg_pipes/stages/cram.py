"""
Stage that generates a CRAM file.
"""

import logging

from .. import Path
from ..filetypes import CramPath
from ..pipeline import stage, SampleStage, StageInput, StageOutput, Sample, PipelineError
from ..jobs import align

logger = logging.getLogger(__file__)


@stage(analysis_type='cram')
class CramStage(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """
    def expected_result(self, sample: Sample) -> Path:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        return sample.get_cram_path().path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "align" function implemented in the jobs module
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

        jobs = align.align(
            b=self.b,
            alignment_input=sample.alignment_input,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            job_attrs=sample.get_job_attrs(),
            refs=self.refs,
            overwrite=not self.check_intermediates,
            number_of_shards_for_realignment=(
                10 if isinstance(sample.alignment_input, CramPath) else None
            ),
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=jobs
        )
