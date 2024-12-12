import logging
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_data_api

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
    """Prepare ICA for data upload and pipeline runs by generating a file ID for
    the CRAM and CRAI files, and a folder ID for the outputs of the Dragen pipeline.

    Args:
        cram (str): The filename of the CRAM that is to be uploaded
        upload_folder (str): The folder in which data should be uploaded into
        ica_analysis_output_folder (str): The folder that outputs from the pipeline run should be written to
        api_root (str): The ICA API endpoint
        sg_name (str): The name of the sequencing group
        bucket_name (str): The  name of the GCP bucket that the data reside in
        gcp_folder (str): The path to the folder in GCP that expected outputs should be written to
    """
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
            output_object_name: str = f'{item["object"]}_ica_file_id'
            if item['object_type'] == 'FOLDER':
                folder_path = f'/{bucket_name}/{ica_analysis_output_folder}'
                output_object_name = f'{item["object"]}_dragen_output_folder_id'
            logging.info(f'File is: {item["object"]}, object type is {item["object_type"]}')
            object_id: str = ica_utils.create_upload_object_id(
                api_instance=api_instance,
                path_params=path_parameters,
                sg_name=sg_name,
                file_name=item['object'],
                folder_path=folder_path,
                object_type=item['object_type'],
            )
            ica_utils.register_output_to_gcp(
                bucket=bucket_name,
                object_contents=object_id,
                object_name=output_object_name,
                gcp_folder=gcp_folder,
            )
