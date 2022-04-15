"""
Stage that generates a CRAM file.
"""

import logging

from .. import Path, types
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from ..jobs import align

logger = logging.getLogger(__file__)


@stage
class Align(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        return sample.get_cram_path().path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "align" function implemented in the jobs module
        """
        if not sample.alignment_input or (
            self.check_inputs and 
            not types.alignment_input_exists(sample.alignment_input)
        ):
            if self.skip_samples_with_missing_input:
                logger.error(f'No alignment inputs, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, 
                    error_msg=f'No alignment input found for {sample.id}'
                )

        jobs = align.align(
            b=self.b,
            alignment_input=sample.alignment_input,
            output_path=self.expected_outputs(sample),
            sample_name=sample.id,
            job_attrs=self.get_job_attrs(sample),
            refs=self.refs,
            overwrite=not self.check_intermediates,
            realignment_shards_num=self.pipeline_config.get('realignment_shards_num'),
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
