import logging
import time
from random import randint
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_analysis_api

from cpg_workflows.stages.dragen_ica import ica_utils


def run(
    ica_pipeline_id_path: str,
    ica_output_folder_id_path: str,
    api_root: str,
    gcp_bucket: str,
    sg_name: str,
    pipeline_registration_path: str,
) -> None:
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key

    ica_pipeline_id: str = ica_utils.read_blob_contents(full_blob_path=ica_pipeline_id_path)
    ica_output_folder_id: str = ica_utils.read_blob_contents(full_blob_path=ica_output_folder_id_path)

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_analysis_api.ProjectAnalysisApi(api_client)
        path_params: dict[str, str] = {'projectId': project_id}

        pipeline_status: str = ica_utils.check_ica_pipeline_status(
            api_instance=api_instance,
            path_params=path_params | {'analysisId': ica_pipeline_id},
        )
        # Other running statuses are REQUESTED AWAITINGINPUT INPROGRESS
        while pipeline_status not in ['SUCCEEDED', 'FAILED', 'FAILEDFINAL', 'ABORTED']:
            time.sleep(600 + randint(-60, 60))
            pipeline_status = ica_utils.check_ica_pipeline_status(
                api_instance=api_instance,
                path_params=path_params | {'analysisId': ica_pipeline_id},
            )
        if pipeline_status == 'SUCCEEDED':
            logging.info(f'Pipeline run {ica_pipeline_id} has succeeded')
            ica_utils.register_output_to_gcp(
                bucket=gcp_bucket,
                object_contents=ica_output_folder_id,
                object_name=f'{sg_name}_pipeline_success',
                gcp_folder=pipeline_registration_path,
            )
        else:
            raise Exception
