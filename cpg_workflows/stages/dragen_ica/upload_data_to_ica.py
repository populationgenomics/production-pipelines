import logging
import subprocess
from typing import Any

import coloredlogs
import icasdk
from google.cloud import storage
from icasdk.apis.tags import project_data_api
from icasdk.model.create_data import CreateData

from cpg_utils.config import get_gcp_project


def check_object_already_exists(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    sg_name: str,
    folder_path: str,
) -> str | None:
    query_params = {
        'filename': [sg_name],
        'filenameMatchMode': 'EXACT',
        'filePath': [f'/{folder_path}/{sg_name}'],
        'filePathMatchMode': 'STARTS_WITH_CASE_INSENSITIVE',
        'type': 'FILE',
    }
    logging.info(f'{query_params}')
    try:
        object_check_response = upload_api_instance.get_project_data_list(
            path_params=path_params,
            query_params=query_params,
        )
        # No data under this filder / sg_name combo, so ok to proceed
        if len(object_check_response.body['items']) == 0:
            return None
        else:
            # TODO We don't check any other statuses, something to be improved in the future.
            # Statuses are ["PARTIAL", "AVAILABLE", "ARCHIVING", "ARCHIVED", "UNARCHIVING", "DELETING", ]
            if object_check_response.body['items'][0]['data']['details']['status'] == 'PARTIAL':
                return object_check_response.body['items'][0]['data']['id']
            raise NotImplementedError('Checking for other status is not implemented yet.')
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> get_project_data_list: {e}') from e


def create_upload_file_id(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    sg_name: str,
    folder_path: str,
) -> str:
    existing_file_id: str | None = None
    logging.info(f'Existing file id before check: {existing_file_id}')
    existing_file_id = check_object_already_exists(upload_api_instance, path_params, sg_name, folder_path)
    logging.info(f'Existing file id after check: {existing_file_id}')
    if not existing_file_id:
        body = CreateData(
            name=sg_name,
            folderPath=f'{folder_path}/',
            dataType='FILE',
        )
        try:
            upload_file_id_response = upload_api_instance.create_data_in_project(
                path_params=path_params,
                body=body,
            )
            return upload_file_id_response.body['data']['id']
        except icasdk.ApiException as e:
            raise icasdk.ApiException(
                f'Exception when calling ProjectDataApi -> create_data_in_project: {e}',
            ) from e
    return existing_file_id


def create_upload_url(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    file_id: str,
) -> str:
    upload_url_path_params: dict[str, str] = path_params | {'dataId': file_id}
    query_params: dict[Any, Any] = {}
    try:
        upload_api_response = upload_api_instance.create_upload_url_for_data(
            path_params=upload_url_path_params,
            query_params=query_params,
        )
        logging.info('Returning URL for upload')
        return upload_api_response.body['url']
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> create_upload_url_for_data: {e}') from e


def upload_data(
    upload_url: str,
    data_to_upload: str,
    bucket: str,
    suffix: str,
) -> None:
    storage_client = storage.Client()

    gcp_bucket = storage_client.bucket(bucket_name=bucket)
    blob_to_download = gcp_bucket.get_blob(f'{suffix}{data_to_upload}')
    blob_to_download.download_to_filename(data_to_upload, timeout=3600)

    logging.info('Uploading data with cURL')
    subprocess.run(['curl', '--upload-file', data_to_upload, f'{upload_url}'])


def register_gcp_output(bucket: str, ica_file_id: str, item: str) -> None:
    upload_name_and_prefix: str = f'ica/{item}'
    blob_client = storage.Blob(name=upload_name_and_prefix, bucket=bucket)
    blob_client.upload_from_string(data=ica_file_id)


def run(
    suffix: str,
    cram: str,
    bucket_name: str,
    upload_folder: str,
    api_root: str,
    project_id: str,
    api_key: str,
):
    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}
    folder_path: str = f'{get_gcp_project()}/{upload_folder}'
    cram_index: str = f'{cram}.crai'

    with icasdk.ApiClient(configuration) as upload_api_client:
        for item in [cram, cram_index]:
            upload_api_instance = project_data_api.ProjectDataApi(upload_api_client)
            logging.info(f'Item is: {item}')
            upload_file_id: str = create_upload_file_id(upload_api_instance, path_parameters, item, folder_path)
            upload_url: str = create_upload_url(upload_api_instance, path_parameters, upload_file_id)
            logging.info(f'Data to upload: {item}')
            upload_data(upload_url, item, bucket_name, suffix)
            register_gcp_output(bucket_name, upload_file_id, item)
