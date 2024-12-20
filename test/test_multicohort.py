"""
Test reading inputs into a Cohort object.
"""

import logging
import re

from pytest_mock import MockFixture

from cpg_workflows.inputs import MultiCohort

from . import set_config

LOGGER = logging.getLogger(__name__)


def _multicohort_config(tmp_path) -> str:
    conf = f"""
    [workflow]
    dataset = 'projecta'
    input_cohorts = ['COH123', 'COH456']

    path_scheme = 'local'

    [storage.default]
    default = '{tmp_path}'

    [storage.projecta]
    default = '{tmp_path}'

    [storage.projectb]
    default = '{tmp_path}'

    [references.broad]
    ref_fasta = 'stub'
    """

    return conf


def mock_get_sgs_by_cohort_project_a(*args, **kwargs) -> list[dict]:
    return [
        {
            'id': 'CPGXXXX',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projecta',
                },
                'externalId': 'NA12340',
                'participant': {
                    'id': 1,
                    'externalId': '8',
                    'reportedSex': 'Male',
                    'meta': {'participant_meta': 'is_here'},
                },
            },
            'assays': [
                {
                    'id': 1,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS28',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGAAAA',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projecta',
                },
                'externalId': 'NA12489',
                'participant': {
                    'id': 2,
                    'externalId': '14',
                    'reportedSex': None,
                    'meta': {'participant_metadata': 'number_fourteen'},
                },
            },
            'assays': [
                {
                    'id': 2,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS29',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def mock_get_sgs_by_cohort_project_b(*args, **kwargs) -> list[dict]:
    return [
        {
            'id': 'CPGCCCCCC',
            'meta': {'sg_meta': 'is_awesome'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projectb',
                },
                'externalId': 'NA111111',
                'participant': {
                    'id': 1,
                    'externalId': '10',
                    'reportedSex': 'Male',
                    'meta': {'participant_meta': 'is_here'},
                },
            },
            'assays': [
                {
                    'id': 1,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS28',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGDDDDDD',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projectb',
                },
                'externalId': 'NA12489',
                'participant': {
                    'id': 2,
                    'externalId': '14',
                    'reportedSex': None,
                    'meta': {'participant_metadata': 'number_fourteen'},
                },
            },
            'assays': [
                {
                    'id': 2,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS29',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def mock_get_sgs_by_cohort_project_a_test_suffix(*args, **kwargs) -> list[dict]:
    return [
        {
            'id': 'CPGXXXX',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projecta-test',
                },
                'externalId': 'NA12340',
                'participant': {
                    'id': 1,
                    'externalId': '8',
                    'reportedSex': 'Male',
                    'meta': {'participant_meta': 'is_here'},
                },
            },
            'assays': [
                {
                    'id': 1,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS28',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGAAAA',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projecta-test',
                },
                'externalId': 'NA12489',
                'participant': {
                    'id': 2,
                    'externalId': '14',
                    'reportedSex': None,
                    'meta': {'participant_metadata': 'number_fourteen'},
                },
            },
            'assays': [
                {
                    'id': 2,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS29',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def mock_get_sgs_by_cohort_project_b_test_suffix(*args, **kwargs) -> list[dict]:
    return [
        {
            'id': 'CPGCCCCCC',
            'meta': {'sg_meta': 'is_awesome'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projectb-test',
                },
                'externalId': 'NA111111',
                'participant': {
                    'id': 1,
                    'externalId': '10',
                    'reportedSex': 'Male',
                    'meta': {'participant_meta': 'is_here'},
                },
            },
            'assays': [
                {
                    'id': 1,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS28',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGDDDDDD',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'technology': 'short-read',
            'type': 'genome',
            'sample': {
                'project': {
                    'name': 'projectb-test',
                },
                'externalId': 'NA12489',
                'participant': {
                    'id': 2,
                    'externalId': '14',
                    'reportedSex': None,
                    'meta': {'participant_metadata': 'number_fourteen'},
                },
            },
            'assays': [
                {
                    'id': 2,
                    'meta': {
                        'platform': '30x Illumina PCR-Free',
                        'concentration': '25',
                        'fluid_x_tube_id': '220405_FS29',
                        'reference_genome': 'Homo sapiens (b37d5)',
                        'volume': '100',
                        'reads_type': 'fastq',
                        'batch': '1',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def mock_get_cohorts(*args, **kwargs) -> dict:
    return {
        'COH123': {
            "sequencing_groups": {
                'projecta': mock_get_sgs_by_cohort_project_a(),
            },
            "name": "CohortA",
        },
        'COH456': {
            "sequencing_groups": {
                'projectb': mock_get_sgs_by_cohort_project_b(),
            },
            "name": "CohortB",
        },
    }


def mock_get_overlapping_cohorts(*args, **kwargs) -> dict:
    return {
        'COH123': {
            "sequencing_groups": {
                'projecta': [mock_get_sgs_by_cohort_project_a()[0]],
            },
            "name": "CohortA",
        },
        'COH456': {
            "sequencing_groups": {
                'projecta': [mock_get_sgs_by_cohort_project_a()[1]],
            },
            "name": "CohortB",
        },
    }


def mock_get_cohorts_with_test_suffix(*args, **kwargs) -> dict:
    return {
        'COH123': {
            "sequencing_groups": {
                'projecta-test': mock_get_sgs_by_cohort_project_a_test_suffix(),
            },
            "name": "CohortA",
        },
        'COH456': {
            "sequencing_groups": {
                'projectb-test': mock_get_sgs_by_cohort_project_b_test_suffix(),
            },
            "name": "CohortB",
        },
    }


def mock_get_analysis_by_sgs(*args, **kwargs) -> dict:
    return {}


def mock_get_pedigree(*args, **kwargs):  # pylint: disable=unused-argument
    return [
        {
            'family_id': 123,
            'individual_id': '8',
            'paternal_id': 14,
            'maternal_id': None,
            'sex': 1,
            'affected': 1,
        },
        {
            'family_id': 124,
            'individual_id': '14',
            'paternal_id': None,
            'maternal_id': None,
            'sex': 2,
            'affected': 1,
        },
    ]


def test_multicohort(mocker: MockFixture, tmp_path):
    """
    Testing creating a Cohort object from metamist mocks.
    """
    set_config(_multicohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)
    mocker.patch('cpg_workflows.metamist.Metamist.get_sgs_for_cohorts', mock_get_cohorts)

    from cpg_workflows.inputs import get_multicohort

    multicohort = get_multicohort()

    assert multicohort
    assert isinstance(multicohort, MultiCohort)

    # Testing Cohort Information
    cohorts = multicohort.get_cohorts()
    assert len(cohorts) == 2
    assert cohorts[0].name == "CohortA"
    assert cohorts[1].name == "CohortB"

    assert len(multicohort.get_sequencing_groups()) == 4
    assert sorted(multicohort.get_sequencing_group_ids()) == ['CPGAAAA', 'CPGCCCCCC', 'CPGDDDDDD', 'CPGXXXX']

    # Test the projects they belong to
    assert multicohort.get_sequencing_groups()[0].dataset.name == 'projecta'
    assert multicohort.get_sequencing_groups()[1].dataset.name == 'projecta'
    assert multicohort.get_sequencing_groups()[2].dataset.name == 'projectb'
    assert multicohort.get_sequencing_groups()[3].dataset.name == 'projectb'

    test_sg_a = multicohort.get_sequencing_groups()[0]
    test_sg_b = multicohort.get_sequencing_groups()[2]

    # Test SequenceGroup Population
    assert test_sg_a.id == 'CPGXXXX'
    assert test_sg_a.external_id == 'NA12340'
    assert test_sg_a.participant_id == '8'
    assert test_sg_a.meta == {'sg_meta': 'is_fun', 'participant_meta': 'is_here', 'phenotypes': {}}

    assert test_sg_b.id == 'CPGCCCCCC'
    assert test_sg_b.external_id == 'NA111111'
    assert test_sg_b.participant_id == '10'
    assert test_sg_b.meta == {'sg_meta': 'is_awesome', 'participant_meta': 'is_here', 'phenotypes': {}}


def test_overlapping_multicohort(mocker: MockFixture, tmp_path):
    """
    Testing multicohorts where different cohorts have sgs from the same dataset.
    """
    set_config(_multicohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)
    mocker.patch('cpg_workflows.metamist.Metamist.get_sgs_for_cohorts', mock_get_overlapping_cohorts)

    from cpg_workflows.inputs import get_multicohort

    multicohort = get_multicohort()

    assert multicohort
    assert isinstance(multicohort, MultiCohort)

    # Testing Cohort Information

    # Validate Cohort Names
    cohorts = multicohort.get_cohorts()
    assert len(cohorts) == 2
    assert cohorts[0].name == "CohortA"
    assert cohorts[1].name == "CohortB"

    # Validate Sequencing Groups
    assert len(multicohort.get_sequencing_groups()) == 2
    assert multicohort.get_sequencing_group_ids() == ['CPGXXXX', 'CPGAAAA']

    # Test the projects they belong to
    assert multicohort.get_sequencing_groups()[0].dataset.name == 'projecta'
    assert multicohort.get_sequencing_groups()[1].dataset.name == 'projecta'

    test_sg_a = multicohort.get_sequencing_groups()[0]

    # Test SequenceGroup Population
    assert test_sg_a.id == 'CPGXXXX'
    assert test_sg_a.external_id == 'NA12340'
    assert test_sg_a.participant_id == '8'
    assert test_sg_a.meta == {'sg_meta': 'is_fun', 'participant_meta': 'is_here', 'phenotypes': {}}

    for cohort in multicohort.get_cohorts():
        assert len(cohort.get_sequencing_group_ids()) == 1


def test_multicohort_dataset_config(mocker: MockFixture, tmp_path):
    """
    This test is to ensure that the dataset name is stripped of the '-test' suffix.
    """
    set_config(_multicohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)
    # mockup the cohort data with '-test' suffix
    mocker.patch('cpg_workflows.metamist.Metamist.get_sgs_for_cohorts', mock_get_cohorts_with_test_suffix)

    from cpg_workflows.inputs import get_multicohort

    multicohort = get_multicohort()

    assert multicohort
    assert isinstance(multicohort, MultiCohort)

    # Testing Cohort Information
    cohorts = multicohort.get_cohorts()
    assert len(cohorts) == 2
    assert cohorts[0].name == "CohortA"
    assert cohorts[1].name == "CohortB"

    assert len(multicohort.get_sequencing_groups()) == 4
    assert sorted(multicohort.get_sequencing_group_ids()) == ['CPGAAAA', 'CPGCCCCCC', 'CPGDDDDDD', 'CPGXXXX']

    # Test the projects they belong to
    # Datasets name should be without '-test' suffix
    assert multicohort.get_sequencing_groups()[0].dataset.name == 'projecta'
    assert multicohort.get_sequencing_groups()[1].dataset.name == 'projecta'
    assert multicohort.get_sequencing_groups()[2].dataset.name == 'projectb'
    assert multicohort.get_sequencing_groups()[3].dataset.name == 'projectb'
