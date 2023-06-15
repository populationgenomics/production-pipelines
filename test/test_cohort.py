"""
Test reading inputs into a Cohort object.
"""

from pytest_mock import MockFixture

from . import set_config


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
                }
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
                }
            ],
        },
    ]


def test_cohort(mocker: MockFixture, tmp_path):
    """
    Testing creating a Cohort object from metamist mocks.
    """

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

    set_config(conf, tmp_path / 'config.toml')

    def mock_get_families(*args, **kwargs):  # pylint: disable=unused-argument
        return []

    def mock_get_pedigree(*args, **kwargs):  # pylint: disable=unused-argument
        return []

    mocker.patch(
        'metamist.apis.FamilyApi.get_families',
        mock_get_families,
    )
    mocker.patch(
        'metamist.apis.FamilyApi.get_pedigree',
        mock_get_pedigree,
    )

    def mock_get_analysis_by_sgs(*args, **kwargs) -> dict:
        return {}

    mocker.patch('cpg_workflows.metamist.Metamist.get_sg_entries', mock_get_sgs)
    mocker.patch(
        'cpg_workflows.metamist.Metamist.get_analyses_by_sid', mock_get_analysis_by_sgs
    )

    # from cpg_workflows.filetypes import BamPath
    from cpg_workflows.inputs import get_cohort

    from cpg_workflows.targets import Sex

    cohort = get_cohort()

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
    assert test_sg.meta == {'sg_meta': 'is_fun', 'participant_meta': 'is_here'}

    # Test Assay Population
    assert test_sg.assays['sequencing'].sequencing_group_id == 'CPGLCL17'
    assert test_sg.assays['sequencing'].id == '1'
    assert test_sg.assays['sequencing'].meta['fluid_x_tube_id'] == '220405_FS28'
    assert (
        test_sg.alignment_input_by_seq_type['genome']
        == test_sg.assays['sequencing'].alignment_input
    )

    assert test_sg.participant_id == '8'
    assert test_sg.pedigree.sex == Sex.MALE
    assert test_sg2.pedigree.sex == Sex.UNKNOWN

    # Test _sequencing_group_by_id attribute
    assert cohort.get_datasets()[0]._sequencing_group_by_id.keys() == {
        'CPGLCL17',
        'CPGLCL25',
    }
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL17'].id == 'CPGLCL17'
    assert cohort.get_datasets()[0]._sequencing_group_by_id['CPGLCL25'].id == 'CPGLCL25'
