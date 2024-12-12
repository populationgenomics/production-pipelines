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
    """Gets the project ID and API key used to interact with ICA

    Raises:
        Exception: Any exception, as we want to fail if we can't get the credentials for any reason

    Returns:
        dict[str, str]: A dictionary with the keys projectId and apiKey
    """
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
    """Check the status of an ICA pipeline via a pipeline ID

    Args:
        api_instance (project_analysis_api.ProjectAnalysisApi): An instance of the ProjectAnalysisApi
        path_params (dict[str, str]): Dict with projectId and analysisId

    Raises:
        icasdk.ApiException: Any exception if the API call is incorrect

    Returns:
        str: The status of the pipeline. Can be one of ['REQUESTED', 'AWAITINGINPUT', 'INPROGRESS', 'SUCCEEDED', 'FAILED', 'FAILEDFINAL', 'ABORTED']
    """
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
    """Check if an object already exists in ICA, as trying to create another object at
    the same path causes an error

    Args:
        api_instance (project_data_api.ProjectDataApi): An instance of the ProjectDataApi
        path_params (dict[str, str]): A dict with the projectId
        file_name (str): The name of the object that you want to check in ICA e.g.
        folder_path (str): The path to the object that you want to create in ICA.
        object_type (str): The type of hte object to create in ICA. Must be one of ['FILE', 'FOLDER']

    Raises:
        NotImplementedError: Only checks for files with the status 'PARTIAL'
        icasdk.ApiException: Other API errors

    Returns:
        str | None: The object ID, if it exists, or else None
    """
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
    """Create an object in ICA that can be used to upload data to,
    or to write analysis outputs into

    Args:
        api_instance (project_data_api.ProjectDataApi): An instance of the ProjectDataApi
        path_params (dict[str, str]): A dict with the projectId
        sg_name (str): The name of the sequencing group
        file_name (str): The name of the file to upload e.g. CPGxxxx.CRAM
        folder_path (str): The base path to the object in ICA to create
        object_type (str): The type of the object to create. Must be one of ['FILE', 'FOLDER']

    Raises:
        icasdk.ApiException: Any API error

    Returns:
        str: The ID of the object that was created, or the existing ID if it was already present.
    """
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
    """Register ICA steps into GCP, so that the workflow can reuse previous steps.

    Args:
        bucket (str): The bucket to write the object in
        object_contents (str): The contents of the object. Usually an ICA ID, but can be any string.
        object_name (str): The name of the object to wwrite.
        gcp_folder (str): The prefix of the object.
    """
    storage_client = storage.Client()
    upload_name_and_prefix: str = f'{gcp_folder}/{object_name}'
    blob_client = storage.Blob(name=upload_name_and_prefix, bucket=storage_client.bucket(bucket_name=bucket))
    blob_client.upload_from_string(data=object_contents)


def read_blob_contents(full_blob_path: str) -> str:
    """Read the contents of a blob in GCP as text

    Args:
        full_blob_path (str): The full path to the blob, including gs://

    Returns:
        str: The contents of the blob as a string. e.g. fil.xxxxx (ICA ID)
    """
    path_components: dict[str, str] = get_path_components_from_gcp_path(full_blob_path)
    gcp_bucket: str = path_components['bucket']
    blob_path: str = f'{path_components["suffix"]}{path_components["file"]}'
    storage_client = storage.Client()
    blob_client = storage.Blob(name=blob_path, bucket=storage_client.bucket(bucket_name=gcp_bucket))
    return blob_client.download_as_text()
