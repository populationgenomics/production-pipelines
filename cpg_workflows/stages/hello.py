"""
Test.
"""

import logging

from cpg_utils import Path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import hello
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@stage
class HelloSequencingGroup(SequencingGroupStage):
    """
    Test.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'bump': seqgroup.dataset.prefix() / 'hello' / f'{seqgroup.id}.bumps',
            'count': seqgroup.dataset.prefix() / 'hello' / f'{seqgroup.id}.counts',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        jobs = hello.per_seqgroup(
            b=get_batch(),
            cram_path=seqgroup.make_cram_path(),
            job_attrs=self.get_job_attrs(seqgroup),
            outputs=self.expected_outputs(seqgroup),
        )
        return self.make_outputs(seqgroup, data=self.expected_outputs(seqgroup), jobs=jobs)


@stage(required_stages=[HelloSequencingGroup])
class HelloCohort(CohortStage):
    """
    Test.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'fred': self.prefix / 'fred.txt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)
        jobs = hello.per_cohort(get_batch(), self.get_job_attrs(cohort), outputs['fred'])
        return self.make_outputs(cohort, data=outputs, jobs=jobs)
