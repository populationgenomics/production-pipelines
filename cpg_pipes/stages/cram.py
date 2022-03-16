"""
Stage that generates a CRAM file.
"""

import logging

from cpg_pipes.jobs import align
from cpg_pipes.pipeline.analysis import AnalysisType, CramPath
from cpg_pipes.pipeline.pipeline import stage, PipelineError
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput

logger = logging.getLogger(__file__)


@stage(sm_analysis_type=AnalysisType.CRAM)
class CramStage(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """
    def expected_result(self, sample: Sample):
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        return sample.get_cram_path().path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "align" function implemented in the jobs module
        """
        if not sample.alignment_input:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(
                    f'No alignment input found for {sample.id}. '
                    f'Checked: Sequence entry and type=CRAM Analysis entry'
                )

        cram_job = align.align(
            b=self.pipe.b,
            alignment_input=sample.alignment_input,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
            overwrite=not self.pipe.check_intermediates,
            smdb=self.pipe.get_db(),
            prev_batch_jobs=self.pipe.prev_batch_jobs,
            number_of_shards_for_realignment=(
                10 if isinstance(sample.alignment_input, CramPath) else None
            )
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[cram_job]
        )
