"""
An example Stage that renames a fastq file to it's external id
"""

from cpg_workflows.workflow import (
    stage,
    Sample,
    SampleStage,
    StageInput,
    StageOutput,
)

from .. import get_batch
from cpg_workflows.jobs import example_job

@stage()
class ExampleStage(SampleStage):
    def expected_outputs(self, sample: Sample):
        # gs://cpg-fewgenomes-test/Workshop22/ExternalId.fastq.gz
        # sample.dataset.prefix() # gs://cpg-fewgenomes-test
        # sample.external_id # ExternalId

        expected_output_path = (
            sample.dataset.prefix() /
            'Workshop22' /
            f'{sample.external_id}.fastq.gz'
        )
        
        return {
            'new_file': expected_output_path
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        expected_output_path = (
            sample.dataset.prefix() /
            'Workshop22' /
            f'{sample.external_id}.fastq.gz'
        )
        input_path = (
            sample.dataset.prefix() / 'WorkshopNov22' / 'BRCA1_R1.fastq.gz'
        )
        job_attrs = {
            'sample': sample.id,
            'samples': [sample.id],
            'external_id': sample.external_id
        }

        jobs = example_job.example_job(
            b=get_batch(),
            output_path=expected_output_path,
            input_path=input_path,
            job_attrs=job_attrs,
        )

        return self.make_outputs(
            sample,
            data=self.expected_outputs(sample),
            jobs=jobs
        )
