"""
Stage that generates a CRAM file.
"""
import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.jobs import align
from cpg_workflows.jobs.align import MissingAlignmentInputException
from cpg_workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)


@stage(analysis_type='cram', analysis_key='cram')
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
            'crai': sample.make_cram_path().index_path,
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
        """
        # We are not listing markduplicates metrics in expected_outputs:
        # if CRAM exist but metrics are for some reason missing, we want to avoid
        # triggering re-running the expensive alignment again. Samtools stats cover
        # the deduplication rate metric. However, we still want to save markduplicates
        # metrics just in case, as they are more detailed.
        markdup_metrics_path = (
            sample.dataset.prefix()
            / 'qc'
            / 'markduplicates_metrics'
            / f'{sample.id}.markduplicates-metrics'
        )

        try:
            jobs = align.align(
                b=get_batch(),
                sample=sample,
                output_path=sample.make_cram_path(),
                job_attrs=self.get_job_attrs(sample),
                overwrite=sample.forced,
                out_markdup_metrics_path=markdup_metrics_path,
            )
            return self.make_outputs(
                sample, data=self.expected_outputs(sample), jobs=jobs
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
