import json
import logging
import time
from random import randint
from typing import Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_analysis_api

import cpg_utils
from cpg_workflows.stages.dragen_ica import ica_utils


def run(
    ica_pipeline_id_path: str,
    api_root: str,
) -> dict[str, str]:
    """Monitor a pipeline running in ICA

    Args:
        ica_pipeline_id_path (str): The path to the JSON file holding the pipeline ID
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

    with open(cpg_utils.to_path(ica_pipeline_id_path), 'rt') as pipeline_fid_handle:
        ica_pipeline_id: str = json.load(pipeline_fid_handle)['pipeline_id']

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
            return {'pipeline': 'success'}
        elif pipeline_status in ['ABORTING', 'ABORTED']:
            logging.info(f'Pipeline run {ica_pipeline_id} has been cancelled')
            return {'pipeline': 'cancelled'}
        elif pipeline_status == 'FAILED':
            return {'pipeline': 'failed'}
        else:
            raise Exception(f'The pipeline run {ica_pipeline_id} has failed, please check ICA for more info.')
