import logging

import icasdk
from icasdk.apis.tags import project_data_api
from icasdk.model.create_data import CreateData

from cpg_workflows.filetypes import CramPath


def create_upload_file_id(
    configuration: icasdk.Configuration,
    path_params: dict[str, str],
    sg_name: str,
    upload_folder: str,
) -> str:
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
