"""
A sample stage.
"""
import logging
from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    Sample,
)

from cpg_workflows.jobs import little_sample_job
from .. import get_batch


@stage()
class TestSampleStage(SampleStage):
    """
    A sample stage.
    """

    def expected_outputs(self, sample: Sample):
        """
        Generate some stuff!
        """
        path = (sample.dataset.prefix() / 'WorkshopNov22' / f'{sample.participant_id}.fastq.gz')
        logging.info(path)
        print(path)
        return {
            'new_file': str(path),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        print(sample.make_cram_path())
        print(self.get_job_attrs(sample))

        logging.info(sample.make_cram_path())

        input_path = (sample.dataset.prefix() / 'Workshop22' / 'BRACA1_R1.fastq.gz')

        jobs = little_sample_job.little_sample_job(
                b=get_batch(),
                output_path=f'gs://cpg-fewgenomes-test/WorkshopNov22/{sample.participant_id}.fastq.gz',
                input_path=input_path,
        )
        return self.make_outputs(
            sample, data=self.expected_outputs(sample), jobs=jobs
        )
