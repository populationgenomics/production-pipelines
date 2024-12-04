import json
import logging
from math import ceil
from typing import TYPE_CHECKING, Final, Literal

import coloredlogs
from google.cloud import secretmanager, storage

import cpg_utils
from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve, get_access_level, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.dragen_ica import run_align_genotype_with_dragen, upload_data_to_ica
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob

GCP_FOLDER_FOR_ICA_UPLOAD: Final = 'ica'
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


@stage(analysis_type='ica_data_upload', analysis_keys=['cram', 'cram_index'])
class UploadDataToIca(SequencingGroupStage):
    # from cpg_workflows.stages.dragen_ica import upload_data_to_ica

    def calculate_needed_storage(
        self,
        cram: str,
        bucket_name: str,
        suffix: str,
    ) -> str:
        storage_client = storage.Client()
        gcp_bucket = storage_client.bucket(bucket_name=bucket_name)
        blob_to_upload_size_bytes: int = gcp_bucket.get_blob(f'{suffix}{cram}').size
        storage_size: int = ceil((blob_to_upload_size_bytes / (1024**3)) + 3)
        return f'{storage_size}Gi'

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, cpg_utils.Path]:
        output_dict: dict[str, cpg_utils.Path] = {
            'cram_id': cpg_utils.to_path(
                f'gs://cpg-{sequencing_group.dataset.name.replace("-test", "")}-{get_access_level()}',
            )
            / GCP_FOLDER_FOR_ICA_UPLOAD
            / f'{sequencing_group.name}.cram_ica_file_id',
            'cram_index_id': cpg_utils.to_path(
                f'gs://cpg-{sequencing_group.dataset.name.replace("-test", "")}-{get_access_level()}',
            )
            / GCP_FOLDER_FOR_ICA_UPLOAD
            / f'{sequencing_group.name}.cram.crai_ica_file_id',
        }
        return output_dict

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        cram_path_components = get_path_components_from_gcp_path(str(sequencing_group.cram))
        suffix: str = cram_path_components['suffix']
        cram: str = cram_path_components['file']
        bucket_name = cram_path_components['bucket']

        upload_job: PythonJob = get_batch().new_python_job(
            name='UploadDataToIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_job.image(image=image_path('cpg_workflows'))

        upload_job.storage(self.calculate_needed_storage(cram, bucket_name, suffix))
        upload_job.call(
            upload_data_to_ica.run,
            suffix=suffix,
            cram=cram,
            bucket_name=bucket_name,
            upload_folder=config_retrieve(['dragen', 'upload_folder']),
            api_root=ICA_REST_ENDPOINT,
            project_id=ICA_PROJECT,
            api_key=API_KEY,
            gcp_folder=GCP_FOLDER_FOR_ICA_UPLOAD,
        )

        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs=upload_job)


@stage(analysis_type='dragen_align_genotype', required_stages=[UploadDataToIca])
class AlignGenotypeWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def read_blob_contents(self, full_blob_path: str) -> str:
        path_components: dict[str, str] = get_path_components_from_gcp_path(full_blob_path)
        gcp_bucket: str = path_components['bucket']
        blob_path: str = f'{path_components["suffix"]}{path_components["file"]}'
        storage_client = storage.Client()
        blob_client = storage.Blob(name=blob_path, bucket=storage_client.bucket(bucket_name=gcp_bucket))
        return blob_client.download_as_text()

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_id: str = self.read_blob_contents(inputs.as_str(sequencing_group, UploadDataToIca, 'cram_id'))
        cram_index_id: str = self.read_blob_contents(inputs.as_str(sequencing_group, UploadDataToIca, 'cram_index_id'))

        dragen_pipeline_id = config_retrieve(['ica.pipelines', 'dragen_3_7_8'])

        dragen_ht_id: str = config_retrieve(['ica.reference_ids', 'dragen_ht_id'])
        cram_reference_id: str = config_retrieve(['ica.reference_ids', 'cram_reference_id'])

        user_tags: list[str] = config_retrieve(['ica.tags', 'user_tags'])
        technical_tags: list[str] = config_retrieve(['ica.tags', 'technical_tags'])
        reference_tags: list[str] = config_retrieve(['ica.tags', 'reference_tags'])
        user_reference: str = config_retrieve(['ica.tags', 'user_reference'])

        align_genotype_job = get_batch().new_python_job(
            name='AlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        align_genotype_job.image(image=image_path('cpg_workflows'))
        align_genotype_job.call(
            run_align_genotype_with_dragen.run,
            cram_id=cram_id,
            cram_index_id=cram_index_id,
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            dragen_pipeline_id=dragen_pipeline_id,
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            api_root=ICA_REST_ENDPOINT,
            project_id=ICA_PROJECT,
            api_key=API_KEY,
        )


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
