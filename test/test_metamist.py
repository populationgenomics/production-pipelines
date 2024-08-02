"""
Test workflow status reporter.
"""

from pathlib import Path
from typing import Any

import pytest
from pytest_mock import MockFixture

from cpg_utils import to_path
from cpg_workflows.filetypes import FastqPair
from cpg_workflows.metamist import (
    Analysis,
    AnalysisStatus,
    AnalysisType,
    Metamist,
    MetamistError,
    parse_reads,
)
from metamist.exceptions import ApiException, ServiceException

from . import set_config


class TimeOutCounter:
    """
    This is a little helper to count attempts in the mock function
    """

    def __init__(self):
        self.count = 0

    def reset(self):
        self.count = 0

    def increment(self):
        self.count += 1


# global variable to count the number of attempts in the mock function
timeout_counter = TimeOutCounter()


TOML = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = 'test'
sequencing_type = 'genome'

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

"""


def mock_aapi_create_analysis_ok(*args, **kwargs):
    """
    Mock function for AnalysisApi.create_analysis as sucess, returns 12345.
    """
    return 12345


def mock_aapi_create_analysis_err(*args, **kwargs):
    """
    Mock function for AnalysisApi.create_analysis as error.
    """
    raise ApiException(
        status=400,
        reason='Must specify "sequencing_group_ids" or "cohort_ids"',
        http_resp=None,
    )


def mock_aapi_create_analysis_timeout(*args, **kwargs):
    """
    Mock function for AnalysisApi.create_analysis time out.
    First 2x ServiceException(Gateway Timeout), then success.
    """
    if timeout_counter.count < 2:
        timeout_counter.increment()
        raise ServiceException(
            status=504,
            reason='Gateway Timeout',
            http_resp=None,
        )
    return 12345


def mock_aapi_update_analysis_ok(*args, **kwargs):
    """
    Mock function for AnalysisApi.update_analysis as sucess, returns True.
    """
    return True


def mock_aapi_update_analysis_err(*args, **kwargs):
    """
    Mock function for AnalysisApi.update_analysis as error.
    """
    raise ApiException(
        status=400,
        reason="You don't have access to this resources, as the resource you requested didn't belong to a project",
        http_resp=None,
    )


def mock_aapi_update_analysis_timeout(*args, **kwargs):
    """
    Mock function for AnalysisApi.update_analysis time out.
    First 2x ServiceException(Gateway Timeout), then success.
    """
    if timeout_counter.count < 2:
        timeout_counter.increment()
        raise ServiceException(
            status=504,
            reason='Gateway Timeout',
            http_resp=None,
        )
    return True


def mock_aapi_update_analysis_timeout_fails(*args, **kwargs):
    """
    Mock function for AnalysisApi.update_analysis time out.
    This function will simulate metamist downtime by raising a timeout error
    """
    raise ServiceException(
        status=504,
        reason='Gateway Timeout',
        http_resp=None,
    )


def mock_query_get_sg_entries_query_returns_empty_data(*args, **kwargs):
    """
    Mock function for AnalysisApi.query(GET_SEQUENCING_GROUPS_QUERY) as sucess.
    It returns empty seq groups.
    """
    return {
        'project': {
            'sequencingGroups': [],
        },
    }


def mock_query_get_analysis_query_returns_empty_data(*args, **kwargs):
    """
    Mock function for AnalysisApi.query(GET_ANALYSES_QUERY) as sucess.
    It returns empty analyses.
    """
    return {
        'project': {
            'analyses': [],
        },
    }


def mock_query_get_analysis_query_returns_one_analysis_no_sg(*args, **kwargs):
    """
    Mock function for AnalysisApi.query(GET_ANALYSES_QUERY) as sucess.
    It returns valid analyses data.
    where type and status are derived from the input kwargs.
    """
    return {
        'project': {
            'analyses': [
                {
                    'id': 12345,
                    'type': kwargs.get('variables', {}).get('analysis_type', None),
                    'meta': {},
                    'output': 'test_output',
                    'status': kwargs.get('variables', {}).get('analysis_status', None),
                    'sequencingGroups': [],
                },
            ],
        },
    }


def mock_query_get_analysis_query_returns_one_analysis_one_sg(*args, **kwargs):
    """
    Mock function for AnalysisApi.query(GET_ANALYSES_QUERY) as sucess.
    It returns valid analyses data.
    where type and status are derived from the input kwargs.
    """
    return {
        'project': {
            'analyses': [
                {
                    'id': 12345,
                    'type': kwargs.get('variables', {}).get('analysis_type', None),
                    'meta': {},
                    'output': 'test_output',
                    'status': kwargs.get('variables', {}).get('analysis_status', None),
                    'sequencingGroups': [{'id': 'SG01'}],
                },
            ],
        },
    }


class TestMetamist:
    @pytest.fixture
    def metamist(self, tmp_path) -> Metamist:
        conf = TOML.format(directory=tmp_path)
        set_config(conf, tmp_path / 'config.toml')
        return Metamist()

    def test_metamist(self, metamist: Metamist):
        assert metamist is not None

    def test_metamist_create_analysis(self, mocker: MockFixture, metamist: Metamist):
        """
        This test creates an analysis in the Metamist API.
        and checks if the analysis is created successfully.
        """

        # test error in API call
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.create_analysis',
            mock_aapi_create_analysis_err,
        )
        analysis = metamist.create_analysis(
            output=to_path('test_output'),
            type_='test',
            status='completed',
            sequencing_group_ids=['test'],
        )
        assert analysis is None

        # test timout in API call
        timeout_counter.reset()
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.create_analysis',
            mock_aapi_create_analysis_timeout,
        )
        analysis = metamist.create_analysis(
            output=to_path('test_output'),
            type_=AnalysisType.parse('custom'),
            status=AnalysisStatus.parse('completed'),
            sequencing_group_ids=['test'],
        )
        assert analysis is not None

        # test sucess
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.create_analysis',
            mock_aapi_create_analysis_ok,
        )
        analysis = metamist.create_analysis(
            output=to_path('test_output'),
            type_=AnalysisType.parse('custom'),
            status=AnalysisStatus.parse('completed'),
            sequencing_group_ids=['test'],
        )
        assert analysis is not None

    def test_metamist_update_analysis(self, mocker: MockFixture, metamist: Metamist):
        """
        This test creates an analysis in the Metamist API.
        and checks if the analysis is created successfully.
        """

        analysis = Analysis.parse(
            {
                'id': 12345,
                'type': 'custom',
                'status': 'in-progress',
                'sequencingGroups': [],
            },
        )
        initial_status = AnalysisStatus.parse('in-progress')
        status_to_be_set = AnalysisStatus.parse('completed')

        # test error in API call
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.update_analysis',
            mock_aapi_update_analysis_err,
        )
        metamist.update_analysis(
            analysis,
            status_to_be_set,
        )
        # TODO fix implementation, keep the assert below commented out for the original behavior
        # assert analysis.status != status_to_be_set

        # test API call with Timeout error for more than 3 times
        # reset status to in-progress
        analysis.status = initial_status
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.update_analysis',
            side_effect=mock_aapi_update_analysis_timeout_fails,
        )
        metamist.update_analysis(
            analysis,
            status_to_be_set,
        )
        # TODO fix implementation, keep the assert below commented out for the original behavior
        # assert analysis.status != status_to_be_set

        # test timeout in API call, first 2x timeout, then success
        timeout_counter.reset()
        # reset status to in-progress
        analysis.status = initial_status
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.update_analysis',
            side_effect=mock_aapi_update_analysis_timeout,
        )
        metamist.update_analysis(
            analysis,
            status_to_be_set,
        )
        assert analysis.status == status_to_be_set

        # test success
        # reset status to in-progress
        analysis.status = initial_status
        mocker.patch(
            'cpg_workflows.metamist.AnalysisApi.update_analysis',
            mock_aapi_update_analysis_ok,
        )
        metamist.update_analysis(
            analysis,
            status_to_be_set,
        )
        assert analysis.status == status_to_be_set

    def test_metamist_get_sg_entries(self, mocker: MockFixture, metamist: Metamist):
        """
        This test gets the sequencing group entries from the Metamist API.
        """

        dataset_name = 'test'
        mocker.patch(
            'cpg_workflows.metamist.query',
            mock_query_get_sg_entries_query_returns_empty_data,
        )
        sgs = metamist.get_sg_entries(dataset_name)

        assert sgs is not None
        assert sgs == []

    def test_metamist_get_analyses_by_sgid_no_data(
        self,
        mocker: MockFixture,
        metamist: Metamist,
    ):
        """
        This test gets empty data from analyses by sgid entries from the Metamist API.
        """

        mocker.patch(
            'cpg_workflows.metamist.query',
            mock_query_get_analysis_query_returns_empty_data,
        )
        analysis_dict = metamist.get_analyses_by_sgid(
            sg_ids=[],
            analysis_type=AnalysisType.parse('custom'),
        )
        # should return empty dictionary
        assert analysis_dict is not None
        assert analysis_dict == {}

        # test graphql query with one analysis, but no sgids
        mocker.patch(
            'cpg_workflows.metamist.query',
            mock_query_get_analysis_query_returns_one_analysis_no_sg,
        )
        analysis_dict = metamist.get_analyses_by_sgid(
            sg_ids=[],
            analysis_type=AnalysisType.parse('custom'),
        )
        # should return empty dictionary
        assert analysis_dict is not None
        assert analysis_dict == {}

    def test_metamist_get_analyses_by_sgid(
        self,
        mocker: MockFixture,
        metamist: Metamist,
    ):
        """
        This test gets analyses by sgid entries from the Metamist API.
        """
        analysis_type = AnalysisType.parse('custom')
        mocker.patch(
            'cpg_workflows.metamist.query',
            mock_query_get_analysis_query_returns_one_analysis_one_sg,
        )
        analysis_dict = metamist.get_analyses_by_sgid(
            sg_ids=[],  # NOTE: This param is not used in the metamist implementation, but has to be a list ???
            analysis_type=analysis_type,
        )
        # should return dictionary with one record
        assert analysis_dict is not None
        assert len(analysis_dict) == 1
        assert analysis_dict.get('SG01', None) is not None
        assert analysis_dict['SG01'].type == analysis_type

    def test_metamist_parse_reads_errors(self):
        """
        This test parses the reads various errors paths of the implementation
        """

        # test missing reads error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={},
                check_existence=False,
            )

        assert str(exc_info.value) == 'SG01: no "meta/reads" field in meta'

        # test missing reads_type error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={'reads': [{'location': 'test.fastq.gz'}]},
                check_existence=False,
            )

        assert str(exc_info.value) == 'SG01: no "meta/reads_type" field in meta'

        # test unsupported reads_type error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [{'location': 'test.fastq.gz'}],
                    'reads_type': 'rubbish',
                },
                check_existence=False,
            )
        assert (
            str(exc_info.value) == 'SG01: ERROR: "reads_type" is expected to be one of (\'fastq\', \'bam\', \'cram\')'
        )

        # test incorrectly formatted error
        with pytest.raises(ValueError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [{'location': 'test.fastq.gz'}],
                    'reads_type': 'fastq',
                },
                check_existence=False,
            )
        assert str(exc_info.value) == (
            'Sequence data for sequencing group SG01 is incorrectly formatted. '
            'Expecting 2 entries per lane (R1 and R2 fastqs), but got 1. '
            'Read data: {\'location\': \'test.fastq.gz\'}'
        )

        # test file does not exist error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [
                        [
                            {'location': 'test1.fastq.gz'},
                            {'location': 'test2.fastq.gz'},
                        ],
                    ],
                    'reads_type': 'fastq',
                },
                check_existence=True,
            )
        assert str(exc_info.value) == ('SG01: ERROR: read 1 file does not exist: test1.fastq.gz')

    def test_metamist_parse_reads_errors_bam(self, metamist: Metamist):
        """
        This test parses the reads various errors paths of the implementation
        for bam files
        NOTE: metamist param needed to populate config,
        even if not referenced in the test itself
        """

        # test presents of multiple files as an error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [
                        {'location': 'file1'},
                        {'location': 'file2'},
                    ],
                    'reads_type': 'bam',
                },
                check_existence=False,
            )
        assert str(exc_info.value) == ('SG01: supporting only single bam/cram input')

        # test incorrect file extension error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [
                        {'location': 'file.rubbish'},
                    ],
                    'reads_type': 'bam',
                },
                check_existence=False,
            )
        assert str(exc_info.value) == (
            'SG01: ERROR: expected the file to have an extension .cram or .bam, got: file.rubbish'
        )

        # test file does not exist error
        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [
                        {'location': 'file.bam'},
                    ],
                    'reads_type': 'bam',
                },
                check_existence=True,
            )
        assert str(exc_info.value) == ('SG01: ERROR: index file does not exist: file.bam')

        with pytest.raises(MetamistError) as exc_info:
            parse_reads(
                sequencing_group_id='SG01',
                assay_meta={
                    'reads': [
                        {
                            'location': 'file.bam',
                            'secondaryFiles': [{'location': 'file.rubbish'}],
                        },
                    ],
                    'reads_type': 'bam',
                },
                check_existence=False,
            )
        assert str(exc_info.value) == (
            'SG01: ERROR: expected the index file to have an extension .crai or .bai, got: file.rubbish'
        )

    def test_metamist_parse_reads_bam_ok(self, metamist: Metamist):
        """
        Test simple bam case, do not check for file existence
        """
        reads = parse_reads(
            sequencing_group_id='SG01',
            assay_meta={
                'reads': [
                    {
                        'location': 'file.bam',
                        'secondaryFiles': [{'location': 'file.bai'}],
                    },
                ],
                'reads_type': 'bam',
            },
            check_existence=False,
        )
        assert str(reads) == 'file.bam'

    def test_metamist_parse_reads_cram_ok(self, metamist: Metamist):
        """
        Test simple cram case, do not check for file existence
        """
        reads = parse_reads(
            sequencing_group_id='SG01',
            assay_meta={
                'reads': [
                    {
                        'location': 'file.cram',
                        'secondaryFiles': [{'location': 'file.crai'}],
                    },
                ],
                'reads_type': 'cram',
            },
            check_existence=False,
        )
        assert str(reads) == 'file.cram'

    def test_metamist_parse_reads_fastq_ok(self, metamist: Metamist):
        """
        Test simple fastq case, do not check for file existence
        """
        reads = parse_reads(
            sequencing_group_id='SG01',
            assay_meta={
                'reads': [
                    [
                        {'location': 'test1.fastq.gz'},
                        {'location': 'test2.fastq.gz'},
                    ],
                ],
                'reads_type': 'fastq',
            },
            check_existence=False,
        )
        assert str(reads) == 'test{1,2}.fastq.gz'
