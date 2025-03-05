import json
import logging
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_analysis_api

import cpg_utils
from cpg_workflows.stages.dragen_ica import ica_utils


def run(ica_pipeline_id_path: str, api_root: str) -> dict[str, str]:
    """Cancel a running ICA pipeline via the API

    Args:
        ica_pipeline_id_path (str): The path to the JSON file holding the pipeline ID
        api_root (str): The root for the ICA API

    Raises:
        icasdk.ApiException: Any API error

    Returns:
        dict[str, str]: A cancelled dict to be recorded in GCP noting that the pipeline was cancelled.
                        Includes a timestamp so that a single cancelled pipeline isn't blocking.
    """
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key

    with open(cpg_utils.to_path(ica_pipeline_id_path), 'rt') as pipeline_fid_handle:
        ica_pipeline_id: str = pipeline_fid_handle.read().rstrip()
    path_parameters: dict[str, str] = {'projectId': project_id} | {'analysisId': ica_pipeline_id}

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_analysis_api.ProjectAnalysisApi(api_client)
        try:
            api_response = api_instance.abort_analysis(
                path_params=path_parameters,
                skip_deserialization=True,
            )
            logging.info(f'Sent cancellation request for ICA analysis: {ica_pipeline_id}')
            return {'cancelled': ica_pipeline_id}
        except icasdk.ApiException as e:
            raise icasdk.ApiException(f'Exception when calling ProjectAnalysisApi->abort_analysis: {e}') from e
