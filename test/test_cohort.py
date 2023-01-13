"""
Test reading inputs into a Cohort object.
"""

import toml
from pytest_mock import MockFixture
from . import results_prefix

TOML = f"""
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
default = '{results_prefix()}'

[storage.fewgenomes]
default = '{results_prefix()}'

[large_cohort]
pop_meta_field = 'Superpopulation name'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'

[references.broad]
ref_fasta = 'stub'
"""


def test_cohort(mocker: MockFixture):
    """
    Testing creating a Cohort object from metamist mocks.
    """
    mocker.patch('cpg_utils.config.get_config', lambda: toml.loads(TOML))

    def mock_get_samples(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> list[dict]:
        return [
            {'id': 'CPG00', 'external_id': 'SAMPLE0'},
            {'id': 'CPG01', 'external_id': 'SAMPLE1'},
            {'id': 'CPG02', 'external_id': 'SAMPLE2'},
            {'id': 'CPG03', 'external_id': 'SAMPLE3'},
            {'id': 'CPG04', 'external_id': 'SAMPLE4'},
        ]

    def get_sequences_by_sample_ids(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> dict:
        return {
            'CPG00': [
                {
                    'id': 0,
                    'sample_id': 'CPG00',
                    'type': 'genome',
                    'status': 'completed',
                    'meta': {'reads': [{'location': 'file.bam'}], 'reads_type': 'bam'},
                },
                {
                    'id': 1,
                    'sample_id': 'CPG00',
                    'type': 'exome',
                    'status': 'completed',
                    'meta': {'reads': [{'location': 'file.bam'}], 'reads_type': 'bam'},
                },
            ],
            'CPG01': [
                {
                    'id': 2,
                    'sample_id': 'CPG01',
                    'type': 'genome',
                    'status': 'completed',
                    'meta': {
                        'reads': [
                            [
                                {'location': 'file.R1.fq.gz'},
                                {'location': 'file.R2.fq.gz'},
                            ]
                        ],
                        'reads_type': 'fastq',
                    },
                }
            ],
            'CPG02': [
                {
                    'id': 3,
                    'sample_id': 'CPG02',
                    'type': 'genome',
                    'status': 'completed',
                    'meta': {'reads': [{'location': 'file.bam'}], 'reads_type': 'bam'},
                }
            ],
            'CPG03': [
                {
                    'id': 4,
                    'sample_id': 'CPG03',
                    'type': 'genome',
                    'status': 'completed',
                    'meta': {'reads': [{'location': 'file.bam'}], 'reads_type': 'bam'},
                }
            ],
            'CPG04': [
                {
                    'id': 5,
                    'sample_id': 'CPG04',
                    'type': 'genome',
                    'status': 'incomplete',
                    'meta': {},
                }
            ],
        }

    def mock_get_external_participant_id_to_internal_sample_id(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> list[list]:
        return [
            ['PART0', 'CPG00'],
            ['PART1', 'CPG01'],
            ['PART2', 'CPG02'],
            ['PART3', 'CPG03'],
            ['PART4', 'CPG04'],
        ]

    def mock_get_participants(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> list[dict]:
        return [
            {
                'external_id': 'PART0',
                'reported_sex': 1,
                'meta': {
                    'Superpopulation name': 'Africa',
                },
            },
            {
                'external_id': 'PART1',
                'reported_sex': 2,
                'meta': {
                    'Dummy': 'dummy',
                },
            },
            {
                'external_id': 'PART2',
                'reported_sex': 0,
            },
            {
                'external_id': 'PART3',
            },
            {
                'external_id': 'PART4',
            },
        ]

    def mock_get_families(*args, **kwargs):  # pylint: disable=unused-argument
        return [{'id': 1}, {'id': 2}]

    def mock_get_pedigree(*args, **kwargs):  # pylint: disable=unused-argument
        return [
            {
                'family_id': 1,
                'individual_id': 'PART1',
                'maternal_id': 'PART2',
                'paternal_id': 0,
                'affected': 2,
                'sex': 2,
            },
            {
                'family_id': 1,
                'individual_id': 'PART2',
                'maternal_id': 0,
                'paternal_id': 0,
                'affected': 1,
                'sex': 2,
            },
            {
                'family_id': 2,
                'individual_id': 'PART3',
                'maternal_id': 0,
                'paternal_id': 0,
                'affected': 0,
                'sex': 0,
            },
        ]

    def mock_query_analyses(*args, **kwargs):  # pylint: disable=unused-argument
        return []

    mocker.patch(
        'sample_metadata.apis.SampleApi.get_samples',
        mock_get_samples,
    )
    mocker.patch(
        'sample_metadata.apis.SequenceApi.get_sequences_by_sample_ids',
        get_sequences_by_sample_ids,
    )
    mocker.patch(
        'sample_metadata.apis.ParticipantApi.get_external_participant_id_to_internal_sample_id',
        mock_get_external_participant_id_to_internal_sample_id,
    )
    mocker.patch(
        'sample_metadata.apis.ParticipantApi.get_participants',
        mock_get_participants,
    )
    mocker.patch(
        'sample_metadata.apis.FamilyApi.get_families',
        mock_get_families,
    )
    mocker.patch(
        'sample_metadata.apis.FamilyApi.get_pedigree',
        mock_get_pedigree,
    )
    mocker.patch(
        'sample_metadata.apis.AnalysisApi.query_analyses',
        mock_query_analyses,
    )

    from cpg_workflows.filetypes import BamPath
    from cpg_workflows.inputs import get_cohort
    from cpg_workflows.targets import Sex

    cohort = get_cohort()
    # the 5th sample doesn't have associated seq/meta/reads
    assert len(cohort.get_samples()) == 5
    assert cohort.get_sample_ids() == ['CPG00', 'CPG01', 'CPG02', 'CPG03', 'CPG04']
    assert cohort.get_samples()[0].id == 'CPG00'
    assert isinstance(
        cohort.get_samples()[0].alignment_input_by_seq_type['genome'], BamPath
    )
    assert isinstance(
        cohort.get_samples()[0].alignment_input_by_seq_type['genome'], BamPath
    )
    assert cohort.get_samples()[0].meta['Superpopulation name'] == 'Africa'
    assert cohort.get_samples()[0].pedigree.sex == Sex.MALE
    assert cohort.get_samples()[0].pedigree.mom is None
    assert cohort.get_samples()[0].pedigree.dad is None
    assert cohort.get_samples()[1].pedigree.sex == Sex.FEMALE
    assert cohort.get_samples()[1].pedigree.mom == cohort.get_samples()[2]
    assert cohort.get_samples()[1].pedigree.dad is None
    assert cohort.get_samples()[1].pedigree.phenotype == 2
    assert cohort.get_samples()[2].pedigree.sex == Sex.FEMALE
    assert cohort.get_samples()[2].pedigree.mom is None
    assert cohort.get_samples()[2].pedigree.dad is None
    assert cohort.get_samples()[3].pedigree.sex == Sex.UNKNOWN
    assert cohort.get_samples()[3].pedigree.mom is None
    assert cohort.get_samples()[3].pedigree.dad is None
    assert cohort.get_samples()[4].seq_by_type['genome'].alignment_input is None
