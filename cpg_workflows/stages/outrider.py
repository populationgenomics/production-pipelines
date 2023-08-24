"""
Perform outlier gene expression analysis with Outrider.
"""

import logging
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
    Cohort,
    CohortStage,
)
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.stages.count import Count
from cpg_workflows.jobs import outrider


@stage(
    required_stages=Count,
)
class Outrider(CohortStage):
    """
    Perform outlier gene expression analysis with Outrider.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate outrider outputs.
        """
        return {
            sequencing_group.id: sequencing_group.dataset.prefix() / 'outrider' / f'{sequencing_group.id}.output'
            for sequencing_group in cohort.get_sequencing_groups()
        }
    
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run outrider.
        """
        count_inputs = [
            inputs.as_path(sequencing_group, Count, 'count')
            for sequencing_group in cohort.get_sequencing_groups()
        ]
        j = outrider.outrider(
            b=get_batch(),
            input_counts=count_inputs,
            output_path=list(self.expected_outputs(cohort).values())[0],
            cohort_name="PLACEHOLDER_NAME",
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)
