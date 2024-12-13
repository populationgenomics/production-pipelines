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


def add_sg(ds, id, external_id: str):
    sg = ds.add_sequencing_group(
        id=id,
        external_id=external_id,
        sequencing_type='genome',
        sequencing_technology='short-read',
        sequencing_platform='illumina',
    )
    return sg


def mock_cohort() -> MultiCohort:
    m = MultiCohort()
    c = m.create_cohort(id='COH123', name='fewgenomes')
    ds = c.create_dataset('my_dataset')
    m_ds = m.add_dataset(ds)

    sg1 = add_sg(ds, 'CPGAA', external_id='SAMPLE1')
    m_ds.add_sequencing_group_object(sg1)
    return m


def mock_multidataset_cohort() -> MultiCohort:
    m = MultiCohort()
    c = m.create_cohort(id='COH123', name='fewgenomes')

    ds = c.create_dataset('my_dataset')
    m_ds_1 = m.add_dataset(ds)

    sg1 = add_sg(ds, 'CPGAA', external_id='SAMPLE1')
    sg2 = add_sg(ds, 'CPGBB', external_id='SAMPLE2')

    m_ds_1.add_sequencing_group_object(sg1)
    m_ds_1.add_sequencing_group_object(sg2)

    ds2 = c.create_dataset('my_dataset2')
    m_ds_2 = m.add_dataset(ds2)

    sg3 = add_sg(ds2, 'CPGCC', external_id='SAMPLE3')
    sg4 = add_sg(ds2, 'CPGDD', external_id='SAMPLE4')
    m_ds_2.add_sequencing_group_object(sg3)
    m_ds_2.add_sequencing_group_object(sg4)

    return m


def mock_multicohort() -> MultiCohort:
    mc = MultiCohort()

    # Create a cohort with two datasets
    cohort_a = mc.create_cohort(id='COH123', name='CohortA')
    # Create a dataset in the cohort (legacy)
    ds = cohort_a.create_dataset('projecta')
    # Create a dataset in the multicohort (new)
    dm1 = mc.add_dataset(ds)
    # Add sequencing groups to the cohort.dataset AND multicohort.dataset
    sg1 = add_sg(ds, 'CPGXXXX', external_id='SAMPLE1')
    sg2 = add_sg(ds, 'CPGAAAA', external_id='SAMPLE2')
    dm1.add_sequencing_group_object(sg1)
    dm1.add_sequencing_group_object(sg2)

    ds2 = cohort_a.create_dataset('projectc')
    dm2 = mc.add_dataset(ds2)
    sg3 = add_sg(ds2, 'CPGCCCC', external_id='SAMPLE3')
    sg4 = add_sg(ds2, 'CPGDDDD', external_id='SAMPLE4')
    dm2.add_sequencing_group_object(sg3)
    dm2.add_sequencing_group_object(sg4)

    cohort_b = mc.create_cohort(id='COH456', name='CohortB')
    ds3 = cohort_b.create_dataset('projectb')
    dm3 = mc.add_dataset(ds3)
    sg5 = add_sg(ds3, 'CPGEEEEEE', external_id='SAMPLE5')
    sg6 = add_sg(ds3, 'CPGFFFFFF', external_id='SAMPLE6')
    dm3.add_sequencing_group_object(sg5)
    dm3.add_sequencing_group_object(sg6)

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
        return to_path(dataset_path(f'{cohort.id}_{self.name}.tsv'))

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
