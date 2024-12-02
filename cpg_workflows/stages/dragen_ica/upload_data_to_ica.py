import logging
from io import BytesIO

import icasdk
import requests
from google.cloud import storage
from icasdk.apis.tags import project_data_api
from icasdk.model.create_data import CreateData

from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import get_gcp_project
from cpg_workflows.filetypes import CramPath


def create_upload_url(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    file_id: str,
) -> str:
    upload_url_path_params: dict[str, str] = path_params | {'dataId': file_id}
    query_params = {}
    try:
        upload_api_response = upload_api_instance.create_upload_url_for_data(
            path_params=upload_url_path_params,
            query_params=query_params,
        )
        return upload_api_response.body['url']
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> create_upload_url_for_data: {e}') from e


def check_object_already_exists(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    sg_name: str,
    folder_path: str,
) -> None:
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
            if object_check_response.body['items'][0]['data']['details']['status'] == 'PARTIAL':
                return object_check_response.body['items'][0]['data']['id']
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> get_project_data_list: {e}') from e


def create_upload_file_id(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    sg_name: str,
    upload_folder: str,
) -> str:
    folder_path: str = f'{get_gcp_project()}/{upload_folder}'
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


def upload_data(upload_url: str, data_to_upload: str, bucket: str) -> None:
    storage_client = storage.Client()

    data_stream = BytesIO()
    gcp_bucket = storage_client.bucket(bucket)
    blob_to_upload = gcp_bucket.blob(data_to_upload)
    blob_to_upload.download_to_file(data_stream)
    blob_to_upload.seek(0)

    res = requests.post(upload_url, data=blob_to_upload)


def run(sg_name: str, sg_path: CramPath, upload_folder: str, api_root: str, project_id: str, api_key: str):
    logging.basicConfig(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}

    # Get the bucket, cram and index information
    cram_path_components = get_path_components_from_gcp_path(str(sg_path.path))
    logging.info(f'{cram_path_components}')
    bucket: str = f'{cram_path_components["bucket"]}'
    cram: str = f'{cram_path_components["suffix"]}{cram_path_components["file"]}'
    cram_index: str = f'{cram_path_components["suffix"]}{cram_path_components["file"]}.crai'

    with icasdk.ApiClient(configuration) as upload_api_client:
        for item in [sg_name, f'{sg_name}.crai']:
            upload_api_instance = project_data_api.ProjectDataApi(upload_api_client)
            upload_file_id: str = create_upload_file_id(upload_api_instance, path_parameters, item, upload_folder)
            upload_url: str = create_upload_url(upload_api_instance, path_parameters, upload_file_id)
            data_to_upload: str = cram if item.endswith('cram') else cram_index

            upload_data(upload_url, data_to_upload, bucket)
