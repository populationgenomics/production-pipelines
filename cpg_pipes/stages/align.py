"""
Stage that generates a CRAM file.
"""

import logging

from cpg_utils.config import get_config

from .. import Path
from ..jobs.align import Aligner, MarkDupTool
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from ..jobs import align
from ..types import CramPath

logger = logging.getLogger(__file__)


@stage(analysis_type='cram')
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
        alignment_input = sample.alignment_input_by_seq_type.get(
            self.cohort.sequencing_type
        )
        if realign_cram_ver := get_config()['workflow'].get('realign_from_cram_version'):
            if (path := (
                sample.dataset.prefix() / 'cram' / realign_cram_ver / f'{sample.id}.cram'
            )).exists():
                logger.info(f'Realigning from {realign_cram_ver} CRAM {path}')
                alignment_input = CramPath(path)

        if alignment_input is None or (
            get_config()['workflow'].get('check_inputs') and not alignment_input.exists()
        ):
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.error(f'No alignment inputs, skipping sample {sample}')
                sample.active = False
                return self.make_outputs(sample, skipped=True)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, error_msg=f'No alignment input found'
                )
        assert alignment_input
        
        jobs = align.align(
            b=self.b,
            alignment_input=alignment_input,
            output_path=self.expected_outputs(sample),
            sample_name=sample.id,
            job_attrs=self.get_job_attrs(sample),
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            realignment_shards_num=get_config()['workflow'].get(
                'realignment_shards_num', align.DEFAULT_REALIGNMENT_SHARD_NUM
            ),
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.PICARD,
            sequencing_type=self.cohort.sequencing_type,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
