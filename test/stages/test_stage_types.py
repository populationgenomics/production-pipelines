"""
Test building stages DAG.
"""

from pytest_mock import MockFixture

from .. import set_config
from . import (
    CohortStage1,
    DatasetStage1,
    DatasetStage2,
    MultiCohortStage1,
    SGStage1,
    TestStage,
    mock_multicohort,
    mock_multidataset_cohort,
    run_workflow,
)


def get_config(tmp_path) -> str:
    conf = f"""
    [workflow]
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    driver_image = 'stub'
    input_datasets = ['fewgenomes']

    check_inputs = false
    check_intermediates = false

    # Skip stages with outputs that already exist
    check_expected_outputs = false

    [storage.default]
    default = '{tmp_path}'
    [storage.fewgenomes]
    default = '{tmp_path}'

    [hail]
    billing_project = 'fewgenomes'
    backend = 'local'
    dry_run = true
    """

    return conf


def get_config_multicohort(tmp_path) -> str:
    conf = f"""
    [workflow]
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'
    driver_image = 'stub'
    input_cohorts = ['COHA', 'COHB']

    check_inputs = false
    check_intermediates = false

    # Skip stages with outputs that already exist
    check_expected_outputs = false

    [storage.default]
    default = '{tmp_path}'
    [storage.fewgenomes]
    default = '{tmp_path}'

    [hail]
    billing_project = 'fewgenomes'
    backend = 'local'
    dry_run = true
    """

    return conf


def test_dataset_stages(mocker: MockFixture, tmp_path):
    """
    Test that dataset stages are run once per dataset.
    """
    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(mocker, cohort_mocker=mock_multidataset_cohort, stages=[DatasetStage1])
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('DatasetStage1')
    # The next two tests are looking at the same behaviour, but from different angles
    assert len(get_multicohort().get_datasets()) == 2
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(get_multicohort().get_datasets())


def test_multiple_dataset_stages(mocker: MockFixture, tmp_path):
    """
    Test that dataset stages are run once per dataset.
    """

    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(mocker, cohort_mocker=mock_multidataset_cohort, stages=[DatasetStage1, DatasetStage2])
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('DatasetStage1')
    assert get_batch().job_by_stage.get('DatasetStage2')

    assert len(get_multicohort().get_datasets()) == 2
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('DatasetStage2').get('job_n') == len(get_multicohort().get_datasets())


def test_mixture_dataset_stage_sg_stage(mocker: MockFixture, tmp_path):
    """
    Tests that when a SequencingGroupStage and DatasetStage are run together, the DatasetStage is run once per dataset and
    the SequencingGroupStage is run once per sequencing group.
    """

    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(mocker, cohort_mocker=mock_multidataset_cohort, stages=[DatasetStage1, SGStage1])
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('DatasetStage1')
    assert get_batch().job_by_stage.get('SGStage1')

    assert len(get_multicohort().get_datasets()) == 2
    assert len(get_multicohort().get_sequencing_groups()) == 4
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('SGStage1').get('job_n') == len(get_multicohort().get_sequencing_groups())


def test_mixture_of_multiple_dataset_and_sg_stages(mocker: MockFixture, tmp_path):
    """
    Tests that when multiple SequencingGroupStages and DatasetStages are run together, the DatasetStage is run once per dataset and
    the SequencingGroupStage is run once per sequencing group.
    """

    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(
        mocker,
        cohort_mocker=mock_multidataset_cohort,
        stages=[
            DatasetStage1,
            DatasetStage2,
            SGStage1,
        ],
    )
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('DatasetStage1')
    assert get_batch().job_by_stage.get('DatasetStage2')
    assert get_batch().job_by_stage.get('SGStage1')

    assert len(get_multicohort().get_datasets()) == 2
    assert len(get_multicohort().get_sequencing_groups()) == 4
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('DatasetStage2').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('SGStage1').get('job_n') == len(get_multicohort().get_sequencing_groups())


def test_depcrecated_cohort_stage_implementation(mocker: MockFixture, tmp_path):
    """
    Tests the CohortStageImplementation, according to the deprecated behaviour (i.e. one cohort per run always)
    """

    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(mocker, cohort_mocker=mock_multidataset_cohort, stages=[CohortStage1])
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('CohortStage1')
    assert get_batch().job_by_stage.get('CohortStage1').get('job_n') == 1

    print(get_multicohort())


def test_mix_all_stage_types(mocker: MockFixture, tmp_path):
    """
    Test CohortStage, DatasetStage, and SequencingGroupStage together.
    """

    set_config(get_config(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort

    run_workflow(
        mocker,
        cohort_mocker=mock_multidataset_cohort,
        stages=[
            CohortStage1,
            DatasetStage1,
            DatasetStage2,
            SGStage1,
        ],
    )
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('CohortStage1')
    assert get_batch().job_by_stage.get('DatasetStage1')
    assert get_batch().job_by_stage.get('DatasetStage2')
    assert get_batch().job_by_stage.get('SGStage1')

    assert len(get_multicohort().get_datasets()) == 2
    assert len(get_multicohort().get_sequencing_groups()) == 4

    assert get_batch().job_by_stage.get('CohortStage1').get('job_n') == 1
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('DatasetStage2').get('job_n') == len(get_multicohort().get_datasets())
    assert get_batch().job_by_stage.get('SGStage1').get('job_n') == len(get_multicohort().get_sequencing_groups())


def test_multicohort_stage(mocker: MockFixture, tmp_path):
    """
    Test MultiCohortStage.
    """

    set_config(get_config_multicohort(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch

    run_workflow(mocker, cohort_mocker=mock_multicohort, stages=[MultiCohortStage1])
    print('Job by stage:', get_batch().job_by_stage)
    assert get_batch().job_by_stage.get('MultiCohortStage1')
    assert get_batch().job_by_stage.get('MultiCohortStage1').get('job_n') == 1


def test_mixture_multicohort_and_other_stages(mocker: MockFixture, tmp_path):
    """
    Testing MultiChortStages with Cohort, Dataset and SequencingGroup stages.
    """

    set_config(get_config_multicohort(tmp_path=tmp_path), tmp_path / 'config.toml')

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import MultiCohort, get_multicohort

    run_workflow(
        mocker,
        cohort_mocker=mock_multicohort,
        stages=[
            CohortStage1,
            DatasetStage1,
            DatasetStage2,
            SGStage1,
            MultiCohortStage1,
        ],
    )
    print('Job by stage:', get_batch().job_by_stage)
    multicohort = get_multicohort()
    assert isinstance(multicohort, MultiCohort)
    assert get_batch().job_by_stage.get('CohortStage1')
    assert get_batch().job_by_stage.get('DatasetStage1')
    assert get_batch().job_by_stage.get('DatasetStage2')
    assert get_batch().job_by_stage.get('SGStage1')
    assert get_batch().job_by_stage.get('MultiCohortStage1')

    assert len(multicohort.get_datasets()) == 3
    assert len(multicohort.get_sequencing_groups()) == 6
    assert len(multicohort.get_cohorts()) == 2

    assert get_batch().job_by_stage.get('MultiCohortStage1').get('job_n') == 1
    assert get_batch().job_by_stage.get('DatasetStage1').get('job_n') == len(multicohort.get_datasets())
    assert get_batch().job_by_stage.get('DatasetStage2').get('job_n') == len(multicohort.get_datasets())
    assert get_batch().job_by_stage.get('SGStage1').get('job_n') == len(multicohort.get_sequencing_groups())
    assert get_batch().job_by_stage.get('CohortStage1').get('job_n') == len(multicohort.get_cohorts())
