"""
Test building stages DAG.
"""
from typing import Callable, Type

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import dataset_path

from cpg_workflows.targets import Sample, Cohort
from cpg_workflows.workflow import (
    SampleStage,
    StageInput,
    StageOutput,
    get_batch,
    run_workflow as _run_workflow,
    stage,
)


def mock_cohort() -> Cohort:
    c = Cohort()

    ds = c.create_dataset('my_dataset')
    ds.add_sample('CPG01', external_id='SAMPLE1')

    return c


class TestStage(SampleStage):
    def expected_outputs(self, sample: Sample) -> Path:
        return to_path(dataset_path(f'{sample.id}_{self.name}.tsv'))

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(sample))
        return self.make_outputs(sample, self.expected_outputs(sample), j)


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
