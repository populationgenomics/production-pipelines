import logging
import subprocess
import time
from random import randint
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_analysis_api

from cpg_workflows.stages.dragen_ica import ica_utils


def run(
    ica_pipeline_id: str | dict[str, str],
    pipeline_id_file: str,
    api_root: str,
) -> dict[str, str]:
    """Monitor a pipeline running in ICA

    Args:
        ica_pipeline_id (str): The path to the file holding the pipeline ID
        api_root (str): The root API endpoint for ICA

    Raises:
        Exception: An exception if the pipeline is cancelled
        Exception: Any other exception if the pipeline gets into a FAILED state

    Returns:
        dict[str, str]: A dict noting success of the pipeline run.
    """
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)

    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key
    pipeline_id: str = ica_pipeline_id['pipeline_id'] if isinstance(ica_pipeline_id, dict) else ica_pipeline_id

    logging.info(f'Monitoring pipeline run {pipeline_id} which of type {type(pipeline_id)}')
    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_analysis_api.ProjectAnalysisApi(api_client)
        path_params: dict[str, str] = {'projectId': project_id}

        pipeline_status: str = ica_utils.check_ica_pipeline_status(
            api_instance=api_instance,
            path_params=path_params | {'analysisId': pipeline_id},
        )
        # Other running statuses are REQUESTED AWAITINGINPUT INPROGRESS
        while pipeline_status not in ['SUCCEEDED', 'FAILED', 'FAILEDFINAL', 'ABORTED']:
            time.sleep(600 + randint(-60, 60))
            pipeline_status = ica_utils.check_ica_pipeline_status(
                api_instance=api_instance,
                path_params=path_params | {'analysisId': pipeline_id},
            )
        if pipeline_status == 'SUCCEEDED':
            logging.info(f'Pipeline run {pipeline_id} has succeeded')
            return {'pipeline': pipeline_id, 'status': 'success'}
        elif pipeline_status in ['ABORTING', 'ABORTED']:
            raise Exception(f'Pipeline run {pipeline_id} has been cancelled.')
        else:
            # Log failed ICA pipeline to a file somewhere
            # Delete the pipeline ID file
            logging.info(f'Deleting the pipeline run ID file {pipeline_id_file}')
            subprocess.run(['gcloud', 'storage', 'rm', pipeline_id_file])
            raise Exception(f'The pipeline run {pipeline_id} has failed, please check ICA for more info.')
