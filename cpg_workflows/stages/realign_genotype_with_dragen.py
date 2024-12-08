import logging
from math import ceil
from typing import TYPE_CHECKING, Final

import coloredlogs
from google.cloud import storage

import cpg_utils
from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve, get_access_level, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.dragen_ica import prepare_ica_for_analysis, run_align_genotype_with_dragen, upload_data_to_ica
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

if TYPE_CHECKING:
    from hailtop.batch.job import PythonJob

GCP_FOLDER_FOR_ICA_UPLOAD: Final = 'ica/prepare'
ICA_REST_ENDPOINT: Final = 'https://ica.illumina.com/ica/rest'


coloredlogs.install(level=logging.INFO)


@stage(analysis_type='prepare_ica_for_analysis', analysis_keys=['cram_fid', 'cram_index_fid', 'analysis_output_fid'])
class PrepareIcaForDragenAnalysis(SequencingGroupStage):
    """Set up ICA for a single realignment run.

    Creates a file ID for both the CRAM and CRAI file to upload to.
    Creates a folder ID for the Dragen output to be written into.

    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, cpg_utils.Path]:
        # This nasty construct is needed in order to stop the pipeline generating a bucket name
        # like fewgenomes-test-test
        sg_bucket: str = f'{sequencing_group.dataset.name.replace("-test", "")}-{get_access_level()}'
        output_dict = {
            'cram_fid': cpg_utils.to_path(
                f'gs://cpg-{sg_bucket}/{GCP_FOLDER_FOR_ICA_UPLOAD}/{sequencing_group.name}.cram_ica_file_id',
            ),
            'cram_index_fid': cpg_utils.to_path(
                f'gs://cpg-{sg_bucket}/{GCP_FOLDER_FOR_ICA_UPLOAD}/{sequencing_group.name}.crai_ica_file_id',
            ),
            'analysis_output_fid': cpg_utils.to_path(
                f'gs://cpg-{sg_bucket}/{GCP_FOLDER_FOR_ICA_UPLOAD}/{sequencing_group.name}.dragen_ouput_folder_id',
            ),
        }
        return output_dict

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path_components = get_path_components_from_gcp_path(str(sequencing_group.cram))
        # suffix: str = cram_path_components['suffix']
        cram: str = cram_path_components['file']
        bucket_name = cram_path_components['bucket']
        logging.info(bucket_name)

        upload_job: PythonJob = get_batch().new_python_job(
            name='UploadDataToIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_job.image(image=image_path('cpg_workflows'))

        upload_job.call(
            prepare_ica_for_analysis.run,
            cram=cram,
            upload_folder=config_retrieve(['dragen', 'upload_folder']),
            ica_analysis_output_folder=config_retrieve(['dragen', 'output_folder']),
            api_root=ICA_REST_ENDPOINT,
            sg_name=sequencing_group.name,
            bucket_name=bucket_name,
            gcp_folder=GCP_FOLDER_FOR_ICA_UPLOAD,
        )

        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs=upload_job)


@stage(
    analysis_type='ica_data_upload',
    analysis_keys=['cram', 'cram_index'],
    required_stages=[PrepareIcaForDragenAnalysis],
)
class UploadDataToIca(SequencingGroupStage):
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
            gcp_folder=GCP_FOLDER_FOR_ICA_UPLOAD,
        )

        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs=upload_job)


@stage(analysis_type='dragen_align_genotype', required_stages=[UploadDataToIca])
class AlignGenotypeWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        return cpg_utils.to_path(f'{sequencing_group.dataset.name}/{sequencing_group.name}/')

    def read_blob_contents(self, full_blob_path: str) -> str:
        path_components: dict[str, str] = get_path_components_from_gcp_path(full_blob_path)
        gcp_bucket: str = path_components['bucket']
        blob_path: str = f'{path_components["suffix"]}{path_components["file"]}'
        storage_client = storage.Client()
        blob_client = storage.Blob(name=blob_path, bucket=storage_client.bucket(bucket_name=gcp_bucket))
        return blob_client.download_as_text()

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_id: str = self.read_blob_contents(str(inputs.as_path(sequencing_group, UploadDataToIca, 'cram_id')))
        cram_index_id: str = self.read_blob_contents(
            str(inputs.as_path(sequencing_group, UploadDataToIca, 'cram_index_id')),
        )

        dragen_pipeline_id = config_retrieve(['ica', 'pipelines', 'dragen_3_7_8'])

        dragen_ht_id: str = config_retrieve(['ica', 'reference_ids', 'dragen_ht_id'])
        cram_reference_id: str = config_retrieve(['ica', 'reference_ids', 'cram_reference_id'])

        user_tags: list[str] = config_retrieve(['ica', 'tags', 'user_tags'])
        technical_tags: list[str] = config_retrieve(['ica', 'tags', 'technical_tags'])
        reference_tags: list[str] = config_retrieve(['ica', 'tags', 'reference_tags'])
        user_reference: str = config_retrieve(['ica', 'tags', 'user_reference'])

        output_folder_path: str = f'{sequencing_group.dataset.name}/{sequencing_group.name}/'

        align_genotype_job = get_batch().new_python_job(
            name='AlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )
        align_genotype_job.image(image=image_path('cpg_workflows'))
        align_genotype_job.call(
            run_align_genotype_with_dragen.run,
            cram_id=cram_id,
            cram_index_id=cram_index_id,
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            dragen_pipeline_id=dragen_pipeline_id,
            output_folder_path=output_folder_path,
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            api_root=ICA_REST_ENDPOINT,
        )
        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs=align_genotype_job)


@stage(analysis_type='dragen_mlr', required_stages=[AlignGenotypeWithDragen])
class GvcfMlrWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        pass


@stage(required_stages=[GvcfMlrWithDragen])
class MonitorGvcfMlrWithDragen(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> None:
        pass

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        pass


@stage(required_stages=[AlignGenotypeWithDragen, GvcfMlrWithDragen])
class CancelIcaPipelineRun(SequencingGroupStage):
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
