import json
import logging
import sys
from typing import TYPE_CHECKING, Final, Literal

import coloredlogs
from google.cloud import secretmanager, storage

from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve, get_access_level, get_cpg_namespace, get_gcp_project, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.dragen_ica import upload_data_to_ica
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob

ICA_REST_ENDPOINT: Final = 'https://ica.illumina.com/ica/rest'
SECRET_CLIENT = secretmanager.SecretManagerServiceClient()
SECRET_PROJECT = 'cpg-common'
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

coloredlogs.install(level=logging.INFO)


@stage(analysis_type='ica_data_upload')
class UploadDataToIca(SequencingGroupStage):
    from cpg_workflows.stages.dragen_ica import upload_data_to_ica

    def expected_outputs(self, sequencing_group: SequencingGroup) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        upload_job: PythonJob = get_batch().new_python_job(
            'UploadDataToIca',
            (self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_job.image(image=image_path('cpg_workflows'))
        # TODO calculate storage for each individual file
        cram_path_components = get_path_components_from_gcp_path(str(sequencing_group.cram))
        cram: str = f'{cram_path_components["suffix"]}{cram_path_components["file"]}'
        bucket_name = cram_path_components['bucket']
        storage_client = storage.Client()
        gcp_bucket = storage_client.bucket(bucket_name=bucket_name)

        logging.info(cram)
        blob_to_upload_size: int | None = gcp_bucket.get_blob(cram).size

        logging.info(blob_to_upload_size)
        sys.exit(0)
        upload_job.storage('40Gi')
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
class AlignGenotypeWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        pass


@stage(analysis_type='dragen_mlr', required_stages=[AlignGenotypeWithDragen])
class GvcfMlrWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        pass


@stage(analysis_type='ica_data_download', required_stages=[AlignGenotypeWithDragen, GvcfMlrWithDragen])
class DownloadDataFromIca(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        pass
