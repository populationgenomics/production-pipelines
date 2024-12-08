import logging
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_data_api

from cpg_utils.config import get_gcp_project
from cpg_workflows.stages.dragen_ica import ica_utils


def run(cram: str, upload_folder: str, ica_analysis_output_folder: str, api_root: str, sg_name: str) -> None:
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}
    folder_path: str = f'{get_gcp_project()}/{upload_folder}'
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
            'object': ica_analysis_output_folder,
            'object_type': 'FOLDER',
        },
    ]
    logging.info('Creating ICA object to upload data.')
    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_data_api.ProjectDataApi(api_client)
        for item in data_setup:
            logging.info(f'File is: {item["object"]}, object type is {item["object_type"]}')
            ica_utils.create_upload_object_id(
                api_instance=api_instance,
                path_params=path_parameters,
                sg_name=sg_name,
                file_name=item['object'],
                folder_path=folder_path,
                object_type=item['object_type'],
            )
