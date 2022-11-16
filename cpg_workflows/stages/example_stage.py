"""
Example sample stage to copy and rename a FASTQ file.
"""
from cpg_workflows.jobs import example_job

from cpg_workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)
from .. import get_batch

@stage
class ExampleStage(SampleStage):
    """
    Copy a FASTQ file and name it using it's external sample ID.
    """

    def expected_outputs(self, sample: Sample):
        """
        Expected output is a FASTQ file named based on the sample's external sample ID.
        """
        expected_output_path = (
            sample.dataset.prefix() / 
            'WorkshopNov22' / 
            f'{sample.external_id}.fastq.gz'
        )

        return {
            'new_file':expected_output_path
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Run the example job function.
        """
        expected_output_path = (
            sample.dataset.prefix() / 
            'WorkshopNov22' / 
            f'{sample.external_id}.fastq.gz'
        )

        input_path = (sample.dataset.prefix() / 'WorkshopNov22' / 'BRCA1_R1.fastq.gz')

        job_attrs = {
            'sample': sample.id,
            'external_id': sample.external_id
        }

        job = example_job.example_job(
            b = get_batch(),
            output_path = expected_output_path,
            input_path = input_path,
            job_attrs = job_attrs,
        )

        return self.make_outputs(sample, data = self.expected_outputs[sample], jobs=[job])
