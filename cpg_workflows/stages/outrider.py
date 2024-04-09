"""
Perform outlier gene expression analysis with Outrider.
"""

from cpg_utils import Path
from cpg_workflows import get_batch
from cpg_workflows.jobs import outrider
from cpg_workflows.stages.count import Count
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(required_stages=Count)
class Outrider(CohortStage):
    """
    Perform outlier gene expression analysis with Outrider.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate outrider outputs.
        """
        dataset_prefix = cohort.get_sequencing_groups()[0].dataset.prefix()
        return {cohort.name: dataset_prefix / 'outrider' / f'{cohort.name}.outrider.RData'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run outrider.
        """
        count_inputs = [
            inputs.as_path(sequencing_group, Count, 'count') for sequencing_group in cohort.get_sequencing_groups()
        ]
        j = outrider.outrider(
            b=get_batch(),
            input_counts=count_inputs,
            output_rdata_path=list(self.expected_outputs(cohort).values())[0],
            cohort_name=cohort.name,
            job_attrs=self.get_job_attrs(),
            overwrite=cohort.forced,
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)
