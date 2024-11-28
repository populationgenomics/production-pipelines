import json
from typing import Final, Literal

from google.cloud import secretmanager

from hailtop.batch.job import PythonJob

from cpg_utils.config import config_retrieve, get_gcp_project, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.dragen_ica import upload_data_to_ica
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import StageOutput, stage

ICA_REST_ENDPOINT: Final = 'https://ica.illumina.com/ica/rest'
SECRET_CLIENT = secretmanager.SecretManagerServiceClient()
SECRET_PROJECT = 'fewgenomes'  # get_gcp_project()  # config_retrieve(['workflow']['project'])  # 'cpg-common'
SECRET_NAME = 'illumina_cpg_workbench_api'
SECRET_VERSION = 'latest'


def get_ica_secrets() -> dict[Literal['projectID', 'apiKey'], str]:
    try:
        secret_path: str = SECRET_CLIENT.secret_version_path(
            project=SECRET_PROJECT,
            secret=SECRET_NAME,
            secret_version=SECRET_VERSION,
        )
        response: secretmanager.AccessSecretVersionResponse = SECRET_CLIENT.access_secret_version(
            request={'name': secret_path},
        )
        return json.loads(response.payload.data.decode('UTF-8'))
    except Exception as e:
        raise Exception(f'Could not obtain ICA credentials: {e}') from e


SECRETS: dict[Literal['projectID', 'apiKey'], str] = get_ica_secrets()
ICA_PROJECT: str = SECRETS['projectID']
API_KEY: str = SECRETS['apiKey']


@stage(anaysis_type='ica_data_upload')
class UploadDataToIca(SequencingGroup):
    from cpg_workflows.stages.dragen_ica import upload_data_to_ica

    def queue_jobs(self, sequencing_group: SequencingGroup) -> StageOutput:
        upload_job: PythonJob = get_batch().new_python_job(
            'UploadDataToIca',
            (self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_job.image(image=image_path('cpg_workflows'))
        upload_job.call(
            upload_data_to_ica.run,
            sg_name=sequencing_group.name,
            sg_path=sequencing_group.cram,
            upload_folder=config_retrieve(['dragen', 'upload_folder']),
            api_root=ICA_REST_ENDPOINT,
            project_id=ICA_PROJECT,
            api_key=API_KEY,
        )


@stage(analysis_type='dragen_align_genotype', required_stages=[UploadDataToIca])
class AlignGenotypeWithDragen(SequencingGroup):
    pass


@stage(analysis_type='dragen_mlr', required_stages=[AlignGenotypeWithDragen])
class GvcfMlrWithDragen(SequencingGroup):
    pass


@stage(analysis_type='ica_data_download')
class DownloadDataFromIca(SequencingGroup):
    pass
