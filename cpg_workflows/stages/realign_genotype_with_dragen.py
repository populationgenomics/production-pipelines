import logging
from math import ceil
from typing import TYPE_CHECKING, Final

import coloredlogs
from google.cloud import storage

from hailtop.batch.job import PythonJob

import cpg_utils
from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import Batch, authenticate_cloud_credentials_in_job, get_batch
from cpg_workflows.stages.dragen_ica import (
    monitor_align_genotype_with_dragen,
    prepare_ica_for_analysis,
    run_align_genotype_with_dragen,
    upload_data_to_ica,
)
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob, PythonJob

DRAGEN_VERSION: Final = config_retrieve(['ica', 'pipelines', 'dragen_version'])
GCP_FOLDER_FOR_ICA_PREP: Final = f'ica/{DRAGEN_VERSION}/prepare'
GCP_FOLDER_FOR_RUNNING_PIPELINE: Final = f'ica/{DRAGEN_VERSION}/pipelines'
GCP_FOLDER_FOR_ICA_DOWNLOAD: Final = f'ica/{DRAGEN_VERSION}/output'
ICA_REST_ENDPOINT: Final = 'https://ica.illumina.com/ica/rest'


coloredlogs.install(level=logging.INFO)


def calculate_needed_storage(
    cram: str,
    bucket_name: str,
) -> str:
    logging.info(f'Checking blob size for {cram}')
    storage_client = storage.Client()
    gcp_bucket = storage_client.bucket(bucket_name=bucket_name)
    blob_to_upload_size_bytes: int = gcp_bucket.get_blob(blob_name=f'{cram}').size
    storage_size: int = ceil(ceil((blob_to_upload_size_bytes / (1024**3)) + 3) * 1.2)
    logging.info(f'Calculated storage is {storage_size}Gi.')
    return f'{storage_size}Gi'


# No need to register this stage in Metamist I think, just ICA prep
@stage(analysis_keys=['cram_fid', 'cram_index_fid', 'analysis_output_fid'])
class PrepareIcaForDragenAnalysis(SequencingGroupStage):
    """Set up ICA for a single realignment run.

    Creates a file ID for both the CRAM and CRAI file to upload to.
    Creates a folder ID for the Dragen output to be written into.

    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, cpg_utils.Path]:
        sg_bucket: str = f'{get_path_components_from_gcp_path(str(sequencing_group.cram))["bucket"]}'
        output_dict = {
            'cram_fid': cpg_utils.to_path(
                f'gs://{sg_bucket}/{GCP_FOLDER_FOR_ICA_PREP}/{sequencing_group.name}.cram_ica_file_id',
            ),
            'crai_fid': cpg_utils.to_path(
                f'gs://{sg_bucket}/{GCP_FOLDER_FOR_ICA_PREP}/{sequencing_group.name}.cram.crai_ica_file_id',
            ),
            'analysis_output_fid': cpg_utils.to_path(
                f'gs://{sg_bucket}/{GCP_FOLDER_FOR_ICA_PREP}/{sequencing_group.name}_dragen_output_folder_id',
            ),
        }
        return output_dict

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path_components: dict[str, str] = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))
        cram: str = cram_path_components['file']
        bucket_name: str = cram_path_components['bucket']
        logging.info(bucket_name)

        prepare_ica_job: PythonJob = get_batch().new_python_job(
            name='PrepareIcaForDragenAnalysis',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        prepare_ica_job.image(image=image_path('cpg_workflows'))

        prepare_ica_job.call(
            prepare_ica_for_analysis.run,
            cram=cram,
            upload_folder=config_retrieve(['ica', 'data_prep', 'upload_folder']),
            ica_analysis_output_folder=config_retrieve(['ica', 'data_prep', 'output_folder']),
            api_root=ICA_REST_ENDPOINT,
            sg_name=sequencing_group.name,
            bucket_name=bucket_name,
            gcp_folder=GCP_FOLDER_FOR_ICA_PREP,
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=prepare_ica_job,
        )


@stage(
    analysis_type='ica_data_upload',
    analysis_keys=['cram_upload_success', 'crai_upload_success'],
    required_stages=[PrepareIcaForDragenAnalysis],
)
class UploadDataToIca(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, cpg_utils.Path]:
        sg_bucket: str = f'{get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))["bucket"]}'
        output_dict: dict[str, cpg_utils.Path] = {
            'cram_upload_success': cpg_utils.to_path(
                f'gs://{sg_bucket}/{GCP_FOLDER_FOR_ICA_PREP}/{sequencing_group.name}.cram_upload_success',
            ),
            'crai_upload_success': cpg_utils.to_path(
                f'gs://{sg_bucket}/{GCP_FOLDER_FOR_ICA_PREP}/{sequencing_group.name}.cram.crai_upload_success',
            ),
        }
        return output_dict

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        cram_path_components: dict[str, str] = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))
        cram: str = f'{cram_path_components["suffix"]}{cram_path_components["file"]}'
        bucket_name: str = cram_path_components['bucket']

        cram_data_mapping: list[dict[str, str]] = [
            {
                'name': f'{cram_path_components["file"]}',
                'full_path': cram,
                'id_path': str(
                    inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis, key='cram_fid'),
                ),
            },
            {
                'name': f'{cram_path_components["file"]}.crai',
                'full_path': f'{cram}.crai',
                'id_path': str(
                    inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis, key='crai_fid'),
                ),
            },
        ]

        upload_job: PythonJob = get_batch().new_python_job(
            name='UploadDataToIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_job.image(image=image_path('cpg_workflows'))

        upload_job.storage(calculate_needed_storage(cram, bucket_name))
        upload_job.call(
            upload_data_to_ica.run,
            cram_data_mapping=cram_data_mapping,
            bucket_name=bucket_name,
            gcp_folder=GCP_FOLDER_FOR_ICA_PREP,
            api_root=ICA_REST_ENDPOINT,
        )

        return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs=upload_job)


@stage(required_stages=[PrepareIcaForDragenAnalysis, UploadDataToIca])
class AlignGenotypeWithDragen(SequencingGroupStage):
    # Output object with pipeline ID to GCP
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        sg_bucket: str = f'{get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))["bucket"]}'
        return cpg_utils.to_path(
            f'gs://{sg_bucket}/{GCP_FOLDER_FOR_RUNNING_PIPELINE}/{sequencing_group.name}_pipeline_id',
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        dragen_pipeline_id = config_retrieve(['ica', 'pipelines', 'dragen_3_7_8'])
        dragen_ht_id: str = config_retrieve(['ica', 'reference_ids', 'dragen_ht_id'])

        cram_reference_id: str = config_retrieve(['ica', 'reference_ids', 'cram_reference_id'])

        user_tags: list[str] = config_retrieve(['ica', 'tags', 'user_tags'])
        technical_tags: list[str] = config_retrieve(['ica', 'tags', 'technical_tags'])
        reference_tags: list[str] = config_retrieve(['ica', 'tags', 'reference_tags'])
        user_reference: str = config_retrieve(['ica', 'tags', 'user_reference'])

        align_genotype_job: PythonJob = get_batch().new_python_job(
            name='AlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )
        align_genotype_job.image(image=image_path('cpg_workflows'))
        align_genotype_job.call(
            run_align_genotype_with_dragen.run,
            cram_name=sequencing_group.name,
            cram_id_path=str(
                inputs.as_path(target=sequencing_group, stage=UploadDataToIca, key='cram_upload_success'),
            ),
            cram_index_id_path=str(
                inputs.as_path(target=sequencing_group, stage=UploadDataToIca, key='crai_upload_success'),
            ),
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            dragen_pipeline_id=dragen_pipeline_id,
            ica_output_folder_id_path=str(
                inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis, key='analysis_output_fid'),
            ),
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            gcp_bucket=get_path_components_from_gcp_path(path=str(sequencing_group.cram))['bucket'],
            pipeline_registration_path=GCP_FOLDER_FOR_RUNNING_PIPELINE,
            api_root=ICA_REST_ENDPOINT,
        )
        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=align_genotype_job,
        )


@stage(analysis_type='dragen_align_genotype', required_stages=[PrepareIcaForDragenAnalysis, AlignGenotypeWithDragen])
class MonitorAlignGenotypeWithDragen(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> cpg_utils.Path:
        sg_bucket: str = f'{get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))["bucket"]}'
        return cpg_utils.to_path(
            f'gs://{sg_bucket}/{GCP_FOLDER_FOR_RUNNING_PIPELINE}/{sequencing_group.name}_pipeline_success',
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        ica_pipeline_id_path: str = str(inputs.as_path(target=sequencing_group, stage=AlignGenotypeWithDragen))
        ica_output_folder_id_path: str = str(
            inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis, key='analysis_output_fid'),
        )

        monitor_pipeline_run: PythonJob = get_batch().new_python_job(
            name='MonitorAlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )

        monitor_pipeline_run.image(image=image_path('cpg_workflows'))
        monitor_pipeline_run.call(
            monitor_align_genotype_with_dragen.run,
            ica_pipeline_id_path=ica_pipeline_id_path,
            ica_output_folder_id_path=ica_output_folder_id_path,
            api_root=ICA_REST_ENDPOINT,
            gcp_bucket=get_path_components_from_gcp_path(path=str(sequencing_group.cram))['bucket'],
            sg_name=sequencing_group.name,
            pipeline_registration_path=GCP_FOLDER_FOR_RUNNING_PIPELINE,
        )
        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=monitor_pipeline_run,
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


@stage(
    analysis_type='ica_data_download',
    required_stages=[PrepareIcaForDragenAnalysis, MonitorAlignGenotypeWithDragen, GvcfMlrWithDragen],
)
class DownloadDataFromIca(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']
        return cpg_utils.to_path(f'gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/{sequencing_group.name}')

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path_components: dict[str, str] = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))
        cram: str = f'{cram_path_components["suffix"]}{cram_path_components["file"]}'
        bucket_name: str = cram_path_components['bucket']
        ica_analysis_folder_id_path: str = str(
            inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis, key='analysis_output_fid'),
        )
        batch_instance: Batch = get_batch()
        ica_download_job: BashJob = batch_instance.new_bash_job(
            name='DownloadDataFromIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        ica_download_job.storage(storage=calculate_needed_storage(cram=cram, bucket_name=bucket_name))
        ica_download_job.image(image=image_path('ica'))

        # Get secrets and folder ID needed for runtime
        ica_analysis_folder_id = batch_instance.read_input(ica_analysis_folder_id_path)

        # Download an entire folder with ICA. Don't log projectId or API key
        authenticate_cloud_credentials_in_job(ica_download_job)
        ica_download_job.command('set +x')
        ica_download_job.command('mkdir -p ~/.icav2')
        ica_download_job.command('echo "server-url: ica.illumina.com" > ~/.icav2/config.yaml')
        ica_download_job.command(
            "gcloud secrets versions access latest --secret=illumina_cpg_workbench_api --project=cpg-common | jq .apiKey | sed 's/\\\"//g' > key",
        )
        ica_download_job.command(
            "gcloud secrets versions access latest --secret=illumina_cpg_workbench_api --project=cpg-common | jq .projectID | sed 's/\\\"//g' > projectID",
        )
        ica_download_job.command(
            'echo "x-api-key: $(cat key)" >> ~/.icav2/config.yaml',
        )
        ica_download_job.command(
            'icav2 projects enter $(cat projectID)',
        )
        ica_download_job.command('set -x')
        ica_download_job.command(
            f'icav2 projectdata download $(cat {ica_analysis_folder_id}) .',
        )
        ica_download_job.command(
            f'mv {bucket_name}/{config_retrieve(["ica", "data_prep", "output_folder"])}/{sequencing_group.name}/* {sequencing_group.name}',
        )
        ica_download_job.command(
            f'cd {sequencing_group.name} && cat *.md5sum > dragen.md5sum && md5sum -c dragen.md5sum',
        )
        ica_download_job.command(
            f'gcloud storage cp --recursive {sequencing_group.name} gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}',
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=ica_download_job,
        )