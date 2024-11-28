import logging

import icasdk
from icasdk.apis.tags import project_data_api
from icasdk.model.create_data import CreateData

from cpg_workflows.filetypes import CramPath


def clean_partial_upload():
    pass


def check_object_already_exists(
    configuration: icasdk.Configuration,
    path_params: dict[str, str],
    sg_name: str,
    upload_folder: str,
) -> None:
    with icasdk.ApiClient(configuration) as object_check_client:
        object_check_instance = project_data_api.ProjectDataApi(object_check_client)
        query_params = {
            'filename': [sg_name],
            'filenameMatchMode': 'EXACT',
            'filePath': [f'/{upload_folder}/{sg_name}/'],
            'filePathMatchMode': 'STARTS_WITH_CASE_INSENSITIVE',
            'type': 'FILE',
        }
        try:
            object_check_response = object_check_instance.get_project_data_list(
                path_params=path_params,
                query_params=query_params,
            )
            # No data under this filder / sg_name combo, so ok to proceed
            if len(object_check_response.body['items']) == 0:
                pass
        except icasdk.ApiException as e:
            raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> get_project_data_list: {e}') from e


def create_upload_file_id(
    configuration: icasdk.Configuration,
    path_params: dict[str, str],
    sg_name: str,
    upload_folder: str,
) -> str:
    check_object_already_exists(configuration, path_params, sg_name, upload_folder)
    with icasdk.ApiClient(configuration) as upload_file_id_client:
        upload_file_id_instance = project_data_api.ProjectDataApi(upload_file_id_client)
        body = CreateData(
            name=sg_name,
            folderPath=upload_folder,
            dataType='FILE',
        )
        try:
            upload_file_id_response = upload_file_id_instance.create_data_in_project(
                path_params=path_params,
                body=body,
            )
            return upload_file_id_response.body['data']['id']
        except icasdk.ApiException as e:
            raise icasdk.ApiException(f'Exception when calling ProjectDataApi -> create_data_in_project: {e}') from e


def run(sg_name: str, sg_path: CramPath, upload_folder: str, api_root: str, project_id: str, api_key: str):
    logging.basicConfig(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}
    upload_file_id: str = create_upload_file_id(configuration, path_parameters, sg_name, upload_folder)
    logging.info(upload_file_id)
