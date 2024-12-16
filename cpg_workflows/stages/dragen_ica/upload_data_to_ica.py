import json
import logging
import subprocess
import sys
from typing import Any, Literal

import coloredlogs
import icasdk
from google.cloud import storage
from icasdk.apis.tags import project_data_api

import cpg_utils
from cpg_workflows.stages.dragen_ica import ica_utils


def create_upload_url(
    upload_api_instance: project_data_api.ProjectDataApi,
    path_params: dict[str, str],
    file_id: str,
) -> str:
    """Generate a presigned URL to upload data to ICA

    Args:
        upload_api_instance (project_data_api.ProjectDataApi): An instance of the ProjectDataApi.
        path_params (dict[str, str]): A Dict of {projectId: id, dataId: id}.
        file_id (str): The ID to populate the path_params dict with.

    Raises:
        icasdk.ApiException: Raises API errors if the API call is formatted incorrectly.

    Returns:
        str: A presigned URL that can be used to upload data.
    """
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
    gcp_path: str,
    object_name: str,
    bucket: str,
) -> None:
    """Uploads data to ICA, via intermediate download to running VM

    Args:
        upload_url (str): The presigned URL to upload the data to
        gcp_path (str): The path in GCP to the object to upload, without gs:// or the bucket name
        object_name (str): The name of the object, used for local download
        bucket (str): The bucket that the object is in in GCP (e.g. fewgenomes-test)
    """
    storage_client = storage.Client()

    gcp_bucket = storage_client.bucket(bucket_name=bucket)
    blob_to_download = gcp_bucket.get_blob(f'{gcp_path}')
    blob_to_download.download_to_filename(object_name, timeout=3600)

    logging.info('Uploading data with cURL')
    subprocess.run(['curl', '--upload-file', object_name, f'{upload_url}'])


def run(
    cram_data_mapping: str,  # list[dict[str, str]],
    bucket_name: str,
    gcp_folder: str,
    api_root: str,
) -> None:
    """Generate a presigned URL per file, and upload the CRAM and CRAI to them.

    Args:
        cram_data_mapping (list[dict[str, str]]): List of dicts, format {name: file_name, full_path: path in gcp minus gs:// and bucket, id_path: GCP path to previous stage output with object ID}
        bucket_name (str): The name of the GCP bucket where the data reside.
        gcp_folder (str): The GCP folder where successful outputs will be written to
        api_root (str): The ICA API endpoint
    """
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']
    coloredlogs.install(level=logging.INFO)

    logging.info(cram_data_mapping)
    logging.info(json.load(open(cpg_utils.to_path(cram_data_mapping))))

    sys.exit(1)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}

    with icasdk.ApiClient(configuration) as upload_api_client:
        upload_api_instance = project_data_api.ProjectDataApi(upload_api_client)
        for item in cram_data_mapping:
            file_id: str = ica_utils.read_blob_contents(full_blob_path=item['id_path'])
            upload_url: str = create_upload_url(
                upload_api_instance=upload_api_instance,
                path_params=path_parameters,
                file_id=file_id,
            )
            upload_data(upload_url=upload_url, gcp_path=item['full_path'], object_name=item['name'], bucket=bucket_name)
            ica_utils.register_output_to_gcp(
                bucket=bucket_name,
                object_contents=file_id,
                object_name=f'{item["name"]}_upload_success',
                gcp_folder=gcp_folder,
            )
