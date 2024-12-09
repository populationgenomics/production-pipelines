import logging
import subprocess
from typing import Any, Literal

import coloredlogs
import icasdk
from google.cloud import storage
from icasdk.apis.tags import project_data_api

from cpg_workflows.stages.dragen_ica import ica_utils


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
) -> None:
    storage_client = storage.Client()

    gcp_bucket = storage_client.bucket(bucket_name=bucket)
    blob_to_download = gcp_bucket.get_blob(f'{data_to_upload}')
    blob_to_download.download_to_filename(data_to_upload, timeout=3600)

    logging.info('Uploading data with cURL')
    subprocess.run(['curl', '--upload-file', data_to_upload, f'{upload_url}'])


def run(
    cram_data_mapping: list[dict[str, str]],
    bucket_name: str,
    gcp_folder: str,
    api_root: str,
):
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']
    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}

    with icasdk.ApiClient(configuration) as upload_api_client:
        upload_api_instance = project_data_api.ProjectDataApi(upload_api_client)
        for item in cram_data_mapping:
            upload_url: str = create_upload_url(upload_api_instance, path_parameters, item['id'])
            upload_data(upload_url, item['file'], bucket_name)
            ica_utils.register_output_to_gcp(bucket_name, 'success', item['name'], gcp_folder)
