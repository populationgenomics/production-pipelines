import logging
from typing import Literal

import coloredlogs


def run(
    ica_analysis_output_folder: str,
    api_root: str,
    sg_name: str,
    bucket_name: str,
) -> dict[str, str]:
    """Prepare ICA pipeline runs by generating a folder ID for the
    outputs of the Dragen pipeline.

    Args:
        ica_analysis_output_folder (str): The folder that outputs from the pipeline run should be written to
        api_root (str): The ICA API endpoint
        sg_name (str): The name of the sequencing group
        bucket_name (str): The  name of the GCP bucket that the data reside in

    Returns:
        dict [str, str] noting the analysis ID.
    """
    import icasdk
    from icasdk.apis.tags import project_data_api

    from cpg_workflows.stages.dragen_ica import ica_utils

    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    path_parameters: dict[str, str] = {'projectId': project_id}

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_data_api.ProjectDataApi(api_client)
        folder_path = f'/{bucket_name}/{ica_analysis_output_folder}'
        object_id: str = ica_utils.create_upload_object_id(
            api_instance=api_instance,
            path_params=path_parameters,
            sg_name=sg_name,
            file_name=sg_name,
            folder_path=folder_path,
            object_type='FOLDER',
        )
        logging.info(f'Created folder ID {object_id} for analysis outputs')
    return {'analysis_output_fid': object_id}
