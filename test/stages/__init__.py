"""
Test building stages DAG.
"""

from typing import Callable, Type

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import dataset_path, get_batch
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)
from cpg_workflows.workflow import (
    run_workflow as _run_workflow,
)


def mock_cohort() -> Cohort:
    c = Cohort()

    ds = c.create_dataset('my_dataset')
    ds.add_sequencing_group('CPG01', external_id='SAMPLE1')

    return c


class TestStage(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> Path:
        return to_path(dataset_path(f'{sequencing_group.id}_{self.name}.tsv'))

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(sequencing_group))
        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), j)


# A -> B -> C -> D
@stage
class A(TestStage):
    pass


@stage(required_stages=A)
class B(TestStage):
    pass


@stage(required_stages=B)
class C(TestStage):
    pass


@stage(required_stages=C)
class D(TestStage):
    pass


# A2 -> B2 -> C2
@stage
class A2(TestStage):
    pass


@stage(required_stages=A2)
class B2(TestStage):
    pass


@stage(required_stages=B2)
class C2(TestStage):
    pass


def run_workflow(
    mocker,
    cohort_mocker: Callable[..., Cohort] = mock_cohort,
    stages: list[Type[TestStage]] | None = None,
):
    mocker.patch('cpg_workflows.inputs.create_cohort', cohort_mocker)

    stages = stages or [C, C2]
    _run_workflow(stages)  # type: ignore
