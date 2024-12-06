import logging
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_data_api

from cpg_utils.config import get_gcp_project
from cpg_workflows.stages.dragen_ica import ica_utils


def run(cram: str, upload_folder: str, ica_analysis_output_folder: str, api_root: str, gcp_folder: str) -> None:
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

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_data_api.ProjectDataApi(api_client)
        for idx, item in enumerate(data_setup):
            for object, object_type in item.items():
                x = ica_utils.create_upload_object_id(
                    api_instance=api_instance,
                    path_params=path_parameters,
                    sg_name=sg_name,
                    file_name=object,
                    folder_path=folder_path,
                    object_type=object_type,
                )
