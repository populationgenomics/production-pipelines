"""
Test workflow status reporter.
"""

from pathlib import Path
from typing import Any

import pytest
from pytest_mock import MockFixture

from cpg_utils import to_path
from cpg_workflows.metamist import Analysis, AnalysisStatus, AnalysisType, Metamist
from metamist.exceptions import ApiException

from . import set_config

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


def mock_query_get_sg_entries_returns_empty_sg_list(*args, **kwargs):
    """
    Mock function for AnalysisApi.get_sg_entries as sucess, returns [].
    """
    return {
        'project': {
            'sequencingGroups': [],
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
        # TODO fix implementation
        # assert analysis.status != status_to_be_set

        # test success
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
            mock_query_get_sg_entries_returns_empty_sg_list,
        )
        sgs = metamist.get_sg_entries(dataset_name)

        assert sgs is not None
        assert sgs == []
