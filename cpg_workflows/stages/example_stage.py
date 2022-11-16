"""
An example stage that doesn't do much 
"""

from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    Sample
)

from cpg_workflows.jobs import example_job
from .. import get_batch 


@stage()
class ExampleStage(SampleStage):
    def expected_outputs(self, sample: Sample):

        expected_output_path = (
            sample.dataset.prefix() /
            'WorkshopNov22' /
            f'{sample.external_id}.fastq.gz'
        )

        return {
            'new_file': expected_output_path
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None :

        expected_output_path = (
            sample.dataset.prefix() / 
            'WorkshopNov22' / 
            f'{sample.external_id}.fastq.gz'
        )

        input_path = (sample.dataset.prefix() / 'WorkshopNov22' / 'BRCA1_R1.fastq.gz')

        job_attrs = {
            'sample': sample.id, 'external_id': sample.external_id
        }

        job = example_job.example_job(
                b=get_batch(),
                output_path=expected_output_path,
                input_path=input_path,
                job_attrs=job_attrs,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[job])
