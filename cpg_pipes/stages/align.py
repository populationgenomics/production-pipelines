"""
Stage that generates a CRAM file.
"""

import logging

from .. import Path
from ..jobs.align import Aligner, MarkDupTool, process_alignment_input
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

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and 
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
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
        assert alignment_input
        
        job_attrs = self.get_job_attrs(sample)
        if seq_type:
            job_attrs['seq_type'] = seq_type.value

        jobs = align.align(
            b=self.b,
            alignment_input=alignment_input,
            output_path=self.expected_outputs(sample),
            sample_name=sample.id,
            job_attrs=self.get_job_attrs(sample),
            overwrite=not self.check_intermediates,
            realignment_shards_num=self.pipeline_config.get(
                'realignment_shards_num', align.DEFAULT_REALIGNMENT_SHARD_NUM
            ),
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.PICARD,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
