"""
Test reading inputs into a Cohort object.
"""

import logging
import re

from pytest_mock import MockFixture

from cpg_workflows.inputs import Cohort, MultiCohort

from . import set_config

LOGGER = logging.getLogger(__name__)


def _cohort_config(tmp_path) -> str:
    conf = f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    sequencing_type = 'genome'

    check_inputs = false
    check_intermediates = false
    check_expected_outputs = false
    path_scheme = 'local'

    [storage.default]
    default = '{tmp_path}'

    [storage.fewgenomes]
    default = '{tmp_path}'

    [large_cohort]
    training_pop = 'Superpopulation name'

    [hail]
    billing_project = 'fewgenomes'
    delete_scratch_on_exit = false
    backend = 'local'

    [references.broad]
    ref_fasta = 'stub'
    """

    return conf


def _custom_cohort_config(tmp_path) -> str:
    conf = f"""
    [workflow]
    dataset = 'projecta'
    input_cohorts = ['COH1']

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


def mock_get_sgs(*args, **kwargs) -> list[dict]:  # pylint: disable=unused-argument
    return [
        {
            'id': 'CPGLCL17',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'reads': [
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1070968,
                                'datetime_added': None,
                            },
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1123158,
                                'datetime_added': None,
                            },
                        ],
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGLCL25',
            'meta': {'sample_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'reads': [
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS29_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS29_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 997128,
                                'datetime_added': None,
                            },
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS29_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS29_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1035385,
                                'datetime_added': None,
                            },
                        ],
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def mock_get_sgs_by_cohort(*args, **kwargs) -> list[dict]:
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


def test_cohort(mocker: MockFixture, tmp_path, caplog):
    """
    Testing creating a Cohort object from metamist mocks.
    """
    set_config(_cohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)

    mocker.patch('cpg_workflows.metamist.Metamist.get_sg_entries', mock_get_sgs)
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)

    caplog.set_level(logging.WARNING)

    from cpg_workflows.inputs import get_multicohort
    from cpg_workflows.targets import SequencingGroup, Sex

    cohort = get_multicohort()

    assert cohort
    assert isinstance(cohort, Cohort)

    # Testing Cohort Information
    assert len(cohort.get_sequencing_groups()) == 2
    assert cohort.get_sequencing_group_ids() == ['CPGLCL17', 'CPGLCL25']

    for sg in cohort.get_sequencing_groups():
        assert sg.dataset.name == 'fewgenomes'
        assert not sg.forced
        assert sg.cram is None
        assert sg.gvcf is None

    # Test SequenceGroup Population
    test_sg = cohort.get_sequencing_groups()[0]
    test_sg2 = cohort.get_sequencing_groups()[1]
    assert test_sg.id == 'CPGLCL17'
    assert test_sg.external_id == 'NA12340'
    assert test_sg.participant_id == '8'
    assert test_sg.meta == {'sg_meta': 'is_fun', 'participant_meta': 'is_here', 'phenotypes': {}}

    # Test Assay Population
    assert test_sg.assays['sequencing'][0].sequencing_group_id == 'CPGLCL17'
    assert test_sg.assays['sequencing'][0].id == '1'
    assert test_sg.assays['sequencing'][0].meta['fluid_x_tube_id'] == '220405_FS28'

    assert test_sg.participant_id == '8'
    # TODO/NOTE: The sex in the pedigree will overwrite the sex in the
    # sequencing group. We should add a check and a test.
    # Also test for unknown reported sex with no ped information.

    assert test_sg.pedigree.sex == Sex.MALE
    assert test_sg2.pedigree.sex == Sex.FEMALE

    assert test_sg.pedigree.mom is None
    assert type(test_sg.pedigree.dad) is SequencingGroup
    assert test_sg.pedigree.dad.participant_id == '14'

    # Test _sequencing_group_by_id attribute
    assert cohort.get_datasets()[0]._sequencing_group_by_id.keys() == {
        'CPGLCL17',
        'CPGLCL25',
    }
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL17'].id == 'CPGLCL17'
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL25'].id == 'CPGLCL25'

    # Test reads
    # TODO: As above
    # assert test_sg.alignment_input_by_seq_type['genome'][0].r1 == CloudPath(
    #     'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz'
    # )
    # assert test_sg2.alignment_input_by_seq_type['genome'][0].r2 == CloudPath(
    #     'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS29_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz'
    # )


def mock_get_sgs_with_missing_reads(*args, **kwargs) -> list[dict]:  # pylint: disable=unused-argument
    return [
        {
            'id': 'CPGLCL17',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'reads': [
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1070968,
                                'datetime_added': None,
                            },
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1123158,
                                'datetime_added': None,
                            },
                        ],
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGLCL25',
            'meta': {'sample_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def test_missing_reads(mocker: MockFixture, tmp_path):
    """
    Testing creating a Cohort object from metamist mocks.
    """
    set_config(_cohort_config(tmp_path), tmp_path / 'config.toml')

    # mock file not existing
    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)

    mocker.patch(
        'cpg_workflows.metamist.Metamist.get_sg_entries',
        mock_get_sgs_with_missing_reads,
    )
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)

    # from cpg_workflows.filetypes import BamPath
    from cpg_workflows.inputs import get_multicohort
    from cpg_workflows.targets import Sex

    cohort = get_multicohort()

    assert cohort

    # Testing Cohort Information
    assert len(cohort.get_sequencing_groups()) == 2
    assert cohort.get_sequencing_group_ids() == ['CPGLCL17', 'CPGLCL25']

    for sg in cohort.get_sequencing_groups():
        assert sg.dataset.name == 'fewgenomes'
        assert not sg.forced
        assert sg.cram is None
        assert sg.gvcf is None

    # Test SequenceGroup Population
    test_sg = cohort.get_sequencing_groups()[0]
    test_sg2 = cohort.get_sequencing_groups()[1]
    assert test_sg.id == 'CPGLCL17'
    assert test_sg.external_id == 'NA12340'
    assert test_sg.participant_id == '8'
    assert test_sg.meta == {'sg_meta': 'is_fun', 'participant_meta': 'is_here', 'phenotypes': {}}

    # Test Assay Population
    assert test_sg.assays['sequencing'][0].sequencing_group_id == 'CPGLCL17'
    assert test_sg.assays['sequencing'][0].id == '1'
    assert test_sg.assays['sequencing'][0].meta['fluid_x_tube_id'] == '220405_FS28'

    assert test_sg.participant_id == '8'
    assert test_sg.pedigree.sex == Sex.MALE
    assert test_sg2.pedigree.sex == Sex.FEMALE

    # Test _sequencing_group_by_id attribute
    assert cohort.get_datasets()[0]._sequencing_group_by_id.keys() == {
        'CPGLCL17',
        'CPGLCL25',
    }
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL17'].id == 'CPGLCL17'
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL25'].id == 'CPGLCL25'

    # Test reads
    # assert test_sg.alignment_input_by_seq_type['genome'][0].r1 == CloudPath(
    #     'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz'
    # )
    assert test_sg2.alignment_input_by_seq_type == {}


def mock_get_sgs_with_mixed_reads(*args, **kwargs) -> list[dict]:  # pylint: disable=unused-argument
    return [
        {
            'id': 'CPGccc',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'reads': [
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1070968,
                                'datetime_added': None,
                            },
                            {
                                'location': 'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'basename': 'HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R2.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1123158,
                                'datetime_added': None,
                            },
                        ],
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGbbb',
            'meta': {'sample_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'genome',
            'sample': {
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
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
        {
            'id': 'CPGaaa',
            'meta': {'sg_meta': 'is_fun'},
            'platform': 'illumina',
            'type': 'exome',
            'sample': {
                'externalId': 'NA1000',
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
                        'reads': [
                            {
                                'location': 'gs://cpg-fewgenomes-main/exomeexample_r1.fastq.gz',
                                'basename': 'exomeexample_r1.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1070968,
                                'datetime_added': None,
                            },
                            {
                                'location': 'gs://cpg-fewgenomes-main/exomeexample_r2.fastq.gz',
                                'basename': 'exomeexample_r2.fastq.gz',
                                'class': 'File',
                                'checksum': None,
                                'size': 1123158,
                                'datetime_added': None,
                            },
                        ],
                        'sequencing_type': 'genome',
                        'sequencing_technology': 'short-read',
                        'sequencing_platform': 'illumina',
                    },
                    'type': 'sequencing',
                },
            ],
        },
    ]


def test_mixed_reads(mocker: MockFixture, tmp_path, caplog):
    """
    Testing creating a Cohort object from metamist mocks.
    """

    caplog.set_level(logging.WARNING)
    set_config(_cohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: True)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)

    mocker.patch(
        'cpg_workflows.metamist.Metamist.get_sg_entries',
        mock_get_sgs_with_mixed_reads,
    )
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)
    from cpg_workflows.inputs import get_multicohort

    cohort = get_multicohort()

    # Testing Cohort Information
    assert len(cohort.get_sequencing_groups()) == 3
    assert cohort.get_sequencing_group_ids() == ['CPGccc', 'CPGbbb', 'CPGaaa']

    # test_genome = cohort.get_sequencing_groups()[0]
    test_none = cohort.get_sequencing_groups()[1]
    # test_exome = cohort.get_sequencing_groups()[2]

    # Test reads
    # TODO: This code returns error: Value of type "AlignmentInput" is not indexable
    # assert test_genome.alignment_input_by_seq_type['genome'][0].r1 == CloudPath(
    #     'gs://cpg-fewgenomes-main/HG3FMDSX3_2_220405_FS28_Homo-sapiens_AACGAGGCCG-ATCCAGGTAT_R_220208_BINKAN1_FEWGENOMES_M001_R1.fastq.gz'
    # )
    # assert test_exome.alignment_input_by_seq_type['exome'][0].r1 == CloudPath(
    #     'gs://cpg-fewgenomes-main/exomeexample_r1.fastq.gz'
    # )

    assert test_none.alignment_input_by_seq_type == {}
    assert re.search(
        r'WARNING\s+root:inputs\.py:\d+\s+No reads found for sequencing group CPGbbb of type genome',
        caplog.text,
    )


def test_unknown_data(mocker: MockFixture, tmp_path, caplog):
    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: True)

    def mock_get_pedigree_empty(*args, **kwargs):
        return []

    def mock_get_families_empty(*args, **kwargs):
        return []

    def mock_get_sgs_with_mixed_reads(*args, **kwargs) -> list[dict]:  # pylint: disable=unused-argument
        return [
            {
                'id': 'CPGccc',
                'meta': {'sg_meta': 'is_fun'},
                'platform': 'illumina',
                'type': 'genome',
                'sample': {
                    'externalId': 'NA12340',
                    'participant': {
                        'id': 1,
                        'externalId': '8',
                        'reportedSex': 'Female',
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
                            'reads': [],
                            'sequencing_type': 'genome',
                            'sequencing_technology': 'short-read',
                            'sequencing_platform': 'illumina',
                        },
                        'type': 'sequencing',
                    },
                ],
            },
            {
                'id': 'CPGaaa',
                'meta': {'sg_meta': 'is_fun'},
                'platform': 'illumina',
                'type': 'exome',
                'sample': {
                    'externalId': 'NA1000',
                    'participant': {
                        'id': 1,
                        'externalId': '8',
                        'reportedSex': None,
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
                            'reads': [
                                {
                                    'location': 'gs://cpg-fewgenomes-main/exomeexample_r1.fastq.gz',
                                    'basename': 'exomeexample_r1.fastq.gz',
                                    'class': 'File',
                                    'checksum': None,
                                    'size': 1070968,
                                    'datetime_added': None,
                                },
                                {
                                    'location': 'gs://cpg-fewgenomes-main/exomeexample_r2.fastq.gz',
                                    'basename': 'exomeexample_r2.fastq.gz',
                                    'class': 'File',
                                    'checksum': None,
                                    'size': 1123158,
                                    'datetime_added': None,
                                },
                            ],
                            'sequencing_type': 'genome',
                            'sequencing_technology': 'short-read',
                            'sequencing_platform': 'illumina',
                        },
                        'type': 'sequencing',
                    },
                ],
            },
        ]

    caplog.set_level(logging.WARNING)
    set_config(_cohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree_empty)

    mocker.patch(
        'cpg_workflows.metamist.Metamist.get_sg_entries',
        mock_get_sgs_with_mixed_reads,
    )
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)
    from cpg_workflows.inputs import get_multicohort
    from cpg_workflows.targets import Sex

    cohort = get_multicohort()

    test_female = cohort.get_sequencing_groups()[0]

    test_unknown = cohort.get_sequencing_groups()[1]

    assert test_female.pedigree.sex == Sex.FEMALE
    assert test_unknown.pedigree.sex == Sex.UNKNOWN


def test_custom_cohort(mocker: MockFixture, tmp_path, monkeypatch):
    """
    Testing creating a Cohort object from metamist mocks.
    """
    set_config(_custom_cohort_config(tmp_path), tmp_path / 'config.toml')

    mocker.patch('cpg_workflows.utils.exists_not_cached', lambda *args: False)

    mocker.patch('cpg_workflows.metamist.Metamist.get_ped_entries', mock_get_pedigree)
    mocker.patch('cpg_workflows.metamist.Metamist.get_analyses_by_sgid', mock_get_analysis_by_sgs)

    def mock_query(query, variables):
        # Mocking the return value of the query function
        return {'cohorts': [{'sequencingGroups': mock_get_sgs_by_cohort()}]}

    # Patching the query function to mock the GraphQL query
    monkeypatch.setattr('cpg_workflows.metamist.query', mock_query)

    from cpg_workflows.inputs import get_multicohort

    cohort = get_multicohort()

    assert cohort
    assert isinstance(cohort, MultiCohort)

    # Testing Cohort Information
    assert len(cohort.get_sequencing_groups()) == 2
    assert cohort.get_sequencing_group_ids() == ['CPGXXXX', 'CPGAAAA']

    # Test the projects they belong to
    assert cohort.get_sequencing_groups()[0].dataset.name == 'projecta'
    assert cohort.get_sequencing_groups()[1].dataset.name == 'projectb'
