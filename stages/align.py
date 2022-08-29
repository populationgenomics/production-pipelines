"""
Stage that generates a CRAM file.
"""
import logging
import dataclasses
from typing import Callable, Optional

from cpg_utils.config import get_config
from cpg_utils import Path
from cpg_utils.workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)

from jobs import align
from jobs.align import Aligner, MarkDupTool, MissingAlignmentInputException


@stage(analysis_type='cram')
class Align(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        return {
            'cram': sample.make_cram_path().path,
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
        """
        jobs = []
        try:
            align_jobs = align.align(
                b=self.b,
                sample=sample,
                output_path=sample.make_cram_path(),
                job_attrs=self.get_job_attrs(sample),
                overwrite=not get_config()['workflow'].get('check_intermediates'),
                aligner=Aligner.DRAGMAP,
                markdup_tool=MarkDupTool.PICARD,
            )
        except MissingAlignmentInputException:
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logging.error(f'No alignment inputs, skipping sample {sample}')
                sample.active = False
                return self.make_outputs(sample, skipped=True)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, error_msg=f'No alignment input found'
                )
        else:
            jobs.extend(align_jobs)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
