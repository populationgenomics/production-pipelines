"""
Test building stages DAG.
"""
from cpg_utils import to_path, Path


def _mock_cohort():
    from cpg_workflows.targets import Cohort

    c = Cohort()
    ds = c.create_dataset('my_dataset')
    ds.add_sample('CPG01', external_id='SAMPLE1')
    return c


def run_workflow(mocker):
    mocker.patch('cpg_workflows.inputs.create_cohort', _mock_cohort)

    from cpg_utils.hail_batch import dataset_path
    from cpg_workflows.targets import Sample
    from cpg_workflows.workflow import (
        SampleStage,
        StageInput,
        StageOutput,
        stage,
        run_workflow,
        get_batch,
    )

    class TestStage(SampleStage):
        def expected_outputs(self, sample: Sample) -> Path:
            return to_path(dataset_path(f'{sample.id}_{self.name}.tsv'))

        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
            j = get_batch().new_job(self.name, attributes=self.get_job_attrs(sample))
            return self.make_outputs(sample, self.expected_outputs(sample), j)

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

    @stage
    class A2(TestStage):
        pass

    @stage(required_stages=A2)
    class B2(TestStage):
        pass

    @stage(required_stages=B2)
    class C2(TestStage):
        pass

    run_workflow([C, C2])
