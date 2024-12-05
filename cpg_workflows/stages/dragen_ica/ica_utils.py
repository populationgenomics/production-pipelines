import json
from typing import Literal

import icasdk
from google.cloud import secretmanager
from icasdk.apis.tags import project_analysis_api


def get_ica_secrets() -> dict[Literal['projectID', 'apiKey'], str]:
    SECRET_CLIENT = secretmanager.SecretManagerServiceClient()
    SECRET_PROJECT = 'cpg-common'
    SECRET_NAME = 'illumina_cpg_workbench_api'
    SECRET_VERSION = 'latest'
    try:
        secret_path: str = SECRET_CLIENT.secret_version_path(
            project=SECRET_PROJECT,
            secret=SECRET_NAME,
            secret_version=SECRET_VERSION,
        )
        response: secretmanager.AccessSecretVersionResponse = SECRET_CLIENT.access_secret_version(
            request={'name': secret_path},
        )
        return json.loads(response.payload.data.decode('UTF-8'))
    except Exception as e:
        raise Exception(f'Could not obtain ICA credentials: {e}') from e


def check_ica_pipeline_status(
    api_instance: project_analysis_api.ProjectAnalysisApi,
    path_params: dict[str, str],
) -> str:
    try:
        api_response = api_instance.get_analysis(path_params=path_params)
        pipeline_status: str = api_response.body['status']
        return pipeline_status
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectAnalysisApi -> get_analysis: {e}') from e
