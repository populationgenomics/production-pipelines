import logging
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_data_api

from cpg_utils.config import get_gcp_project
from cpg_workflows.stages.dragen_ica import ica_utils


def run(
    cram: str,
    upload_folder: str,
    ica_analysis_output_folder: str,
    api_root: str,
    sg_name: str,
    bucket_name: str,
    gcp_folder: str,
) -> None:
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}
    cram_index: str = f'{cram}.crai'
    data_setup: list[dict[str, str]] = [
        {
            'object': cram,
            'object_type': 'FILE',
        },
        {
            'object': cram_index,
            'object_type': 'FILE',
        },
        {
            'object': sg_name,
            'object_type': 'FOLDER',
        },
    ]
    logging.info('Creating ICA object to upload data.')
    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_data_api.ProjectDataApi(api_client)
        for item in data_setup:
            folder_path: str = f'/{bucket_name}/{upload_folder}'
            if item['object_type'] == 'FOLDER':
                folder_path = f'/{bucket_name}/{ica_analysis_output_folder}'
            logging.info(f'File is: {item["object"]}, object type is {item["object_type"]}')
            object_id: str = ica_utils.create_upload_object_id(
                api_instance=api_instance,
                path_params=path_parameters,
                sg_name=sg_name,
                file_name=item['object'],
                folder_path=folder_path,
                object_type=item['object_type'],
            )
            ica_utils.register_output_to_gcp(bucket_name, object_id, f'{item["object"]}_ica_file_id', gcp_folder)
