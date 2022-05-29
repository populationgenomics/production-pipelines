"""
Stage that generates a CRAM file.
"""

import logging

from .. import Path, types
from ..jobs.align import Aligner, MarkDupTool
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from ..jobs import align
from ..types import CramPath, AlignmentInput

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
        seq_type = self.pipeline_config.get('sequencing_type')
        avail_seq_types = set(
            st.value for st in sample.alignment_input_by_seq_type.keys()
        )
        if len(avail_seq_types) > 1 and not seq_type:
            raise ValueError(
                f'{sample}: found alignment inputs with more than one sequencing '
                f'type: {", ".join(avail_seq_types)}. Consider option '
                f'--sequencing-type to use data of specific sequencing type.'
            )
        if seq_type:
            alignment_input = sample.alignment_input_by_seq_type.get(seq_type)
        else:
            alignment_input = list(sample.alignment_input_by_seq_type.values())[0]
            seq_type = alignment_input.sequencing_type

        # Check CRAM for realignment:
        if cram_ver := self.pipeline_config.get('realign_from_cram_version'):
            older_cram = sample.dataset.path() / 'cram' / cram_ver / f'{sample.id}.cram'
            if older_cram.exists():
                logger.info(f'Realigning from {cram_ver} CRAM {older_cram}')
                alignment_input = AlignmentInput(CramPath(older_cram), seq_type)

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

        jobs = align.align(
            b=self.b,
            alignment_input=alignment_input,
            output_path=self.expected_outputs(sample),
            sample_name=sample.id,
            job_attrs=self.get_job_attrs(sample),
            refs=self.refs,
            images=self.images,
            overwrite=not self.check_intermediates,
            realignment_shards_num=self.pipeline_config.get('realignment_shards_num'),
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.PICARD,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
