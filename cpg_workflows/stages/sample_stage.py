"""
A sample stage.
"""

from cpg_utils import to_path
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
        Generate some stuff!.
        """
        return {
            'new_file': f'{sample.participant_id}.cram',
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        print(sample.make_cram_path())
        print(self.get_job_attrs(sample))

        jobs = little_sample_job.little_sample_job(
                b=get_batch(),
                output_path=f'{sample.participant_id}.cram',
                input_path=f'{sample.id}.cram',
        )
        return self.make_outputs(
            sample, data=self.expected_outputs(sample), jobs=jobs
        )
