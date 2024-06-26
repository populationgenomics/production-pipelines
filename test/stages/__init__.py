"""
Test building stages DAG.
"""

from typing import Callable, Type, Union

from cpg_utils import Path, to_path
from cpg_utils.config import dataset_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import Cohort, Dataset, MultiCohort, SequencingGroup
from cpg_workflows.workflow import (
    CohortStage,
    DatasetStage,
    MultiCohortStage,
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
    ds.add_sequencing_group('CPGAA', external_id='SAMPLE1')

    return c


def mock_multidataset_cohort() -> Cohort:
    c = Cohort()

    ds = c.create_dataset('my_dataset')
    ds.add_sequencing_group('CPGAA', external_id='SAMPLE1')
    ds.add_sequencing_group('CPGBB', external_id='SAMPLE2')

    ds2 = c.create_dataset('my_dataset2')
    ds2.add_sequencing_group('CPGCC', external_id='SAMPLE3')
    ds2.add_sequencing_group('CPGDD', external_id='SAMPLE4')

    return c


def mock_multicohort() -> MultiCohort:
    mc = MultiCohort()

    c = mc.create_cohort('CohortA')
    ds = c.create_dataset('projecta')
    ds.add_sequencing_group('CPGXXXX', external_id='SAMPLE1')
    ds.add_sequencing_group('CPGAAAA', external_id='SAMPLE2')

    ds2 = c.create_dataset('projectc')
    ds2.add_sequencing_group('CPGCCCC', external_id='SAMPLE3')
    ds2.add_sequencing_group('CPGDDDD', external_id='SAMPLE4')

    d = mc.create_cohort('CohortB')
    ds3 = d.create_dataset('projectb')
    ds3.add_sequencing_group('CPGEEEEEE', external_id='SAMPLE5')
    ds3.add_sequencing_group('CPGFFFFFF', external_id='SAMPLE6')

    return mc


class TestStage(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> Path:
        return to_path(dataset_path(f'{sequencing_group.id}_{self.name}.tsv'))

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(sequencing_group))
        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), j)


class TestDatasetStage(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> Path:
        return to_path(dataset_path(f'{dataset.name}_{self.name}.tsv'))

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(dataset))
        return self.make_outputs(dataset, self.expected_outputs(dataset), j)


class TestCohortStage(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return to_path(dataset_path(f'{cohort.name}_{self.name}.tsv'))

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(cohort))
        return self.make_outputs(cohort, self.expected_outputs(cohort), j)


class TestMultiCohortStage(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return to_path(dataset_path(f'{multicohort.name}_{self.name}.tsv'))

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        j = get_batch().new_job(self.name, attributes=self.get_job_attrs(multicohort))
        return self.make_outputs(multicohort, self.expected_outputs(multicohort), j)


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


@stage()
class SGStage1(TestStage):
    pass


@stage()
class DatasetStage1(TestDatasetStage):
    pass


@stage()
class DatasetStage2(TestDatasetStage):
    pass


@stage()
class CohortStage1(TestCohortStage):
    pass


@stage()
class MultiCohortStage1(TestMultiCohortStage):
    pass


StageType = Union[Type[TestStage], Type[TestDatasetStage], Type[TestCohortStage], Type[TestMultiCohortStage]]


def run_workflow(
    mocker,
    cohort_mocker: Callable[..., Cohort | MultiCohort] = mock_cohort,
    stages: list[StageType] | None = None,
):
    mocker.patch('cpg_workflows.inputs.deprecated_create_cohort', cohort_mocker)
    mocker.patch('cpg_workflows.inputs.actual_get_multicohort', cohort_mocker)

    stages = stages or [C, C2]
    _run_workflow(stages)  # type: ignore
