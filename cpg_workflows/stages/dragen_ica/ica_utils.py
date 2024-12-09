import json
import logging
from typing import TYPE_CHECKING, Literal

import coloredlogs
import icasdk
from google.cloud import secretmanager, storage
from icasdk.apis.tags import project_analysis_api, project_data_api
from icasdk.model.create_data import CreateData

from cpg_utils.cloud import get_path_components_from_gcp_path

if TYPE_CHECKING:
    from collections.abc import Sequence


coloredlogs.install(level=logging.INFO)


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


def check_object_already_exists(
    api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    file_name: str,
    folder_path: str,
    object_type: str,
) -> str | None:
    query_params: dict[str, Sequence[str] | list[str] | str] = {
        'filePath': [f'{folder_path}/{file_name}'],
        'filePathMatchMode': 'STARTS_WITH_CASE_INSENSITIVE',
        'type': object_type,
    }
    if object_type == 'FILE':
        query_params = {
            'filename': [file_name],
            'filenameMatchMode': 'EXACT',
        } | query_params
    logging.info(f'{query_params}')
    logging.info(f'Checking to see if the {object_type} object already exists at {folder_path}/{file_name}')
    try:
        api_response = api_instance.get_project_data_list(
            path_params=path_params,
            query_params=query_params,
        )
        if len(api_response.body['items']) == 0:
            return None
        elif object_type == 'FOLDER' or api_response.body['items'][0]['data']['details']['status'] == 'PARTIAL':
            return api_response.body['items'][0]['data']['id']
        # Statuses are ["PARTIAL", "AVAILABLE", "ARCHIVING", "ARCHIVED", "UNARCHIVING", "DELETING", ]
        raise NotImplementedError('Checking for other status is not implemented yet.')
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> get_project_data_list: {e}') from e


def create_upload_object_id(
    api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    sg_name: str,
    file_name: str,
    folder_path: str,
    object_type: str,
) -> str:
    existing_object_id: str | None = check_object_already_exists(
        api_instance=api_instance,
        path_params=path_params,
        file_name=file_name,
        folder_path=folder_path,
        object_type=object_type,
    )
    logging.info(f'{existing_object_id}')
    if existing_object_id:
        return existing_object_id
    try:
        if object_type == 'FILE':
            body = CreateData(
                name=file_name,
                folderPath=f'{folder_path}/',
                dataType=object_type,
            )
        else:
            body = CreateData(
                name=sg_name,
                folderPath=f'{folder_path}/',
                dataType=object_type,
            )
        api_response = api_instance.create_data_in_project(
            path_params=path_params,
            body=body,
        )
        return api_response.body['data']['id']
    except icasdk.ApiException as e:
        raise icasdk.ApiException(
            f'Exception when calling ProjectDataApi -> create_data_in_project: {e}',
        ) from e


def register_output_to_gcp(bucket: str, object_contents: str, object_name: str, gcp_folder: str) -> None:
    storage_client = storage.Client()
    upload_name_and_prefix: str = f'{gcp_folder}/{object_name}'
    blob_client = storage.Blob(name=upload_name_and_prefix, bucket=storage_client.bucket(bucket_name=bucket))
    blob_client.upload_from_string(data=object_contents)


def read_blob_contents(full_blob_path: str) -> str:
    path_components: dict[str, str] = get_path_components_from_gcp_path(full_blob_path)
    gcp_bucket: str = path_components['bucket']
    blob_path: str = f'{path_components["suffix"]}{path_components["file"]}'
    storage_client = storage.Client()
    blob_client = storage.Blob(name=blob_path, bucket=storage_client.bucket(bucket_name=gcp_bucket))
    return blob_client.download_as_text()
