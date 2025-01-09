import logging
from datetime import datetime
from math import ceil
from typing import TYPE_CHECKING, Final

import coloredlogs

from hailtop.batch.job import PythonJob

import cpg_utils
from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import Batch, authenticate_cloud_credentials_in_job, get_batch
from cpg_workflows.stages.dragen_ica import (
    cancel_ica_pipeline_run,
    monitor_align_genotype_with_dragen,
    prepare_ica_for_analysis,
    run_align_genotype_with_dragen,
)
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import slugify
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob, PythonJob

DRAGEN_VERSION: Final = config_retrieve(['ica', 'pipelines', 'dragen_version'])
GCP_FOLDER_FOR_ICA_PREP: Final = f'ica/{DRAGEN_VERSION}/prepare'
GCP_FOLDER_FOR_RUNNING_PIPELINE: Final = f'ica/{DRAGEN_VERSION}/pipelines'
GCP_FOLDER_FOR_ICA_DOWNLOAD: Final = f'ica/{DRAGEN_VERSION}/output'
ICA_REST_ENDPOINT: Final = 'https://ica.illumina.com/ica/rest'
ICA_CLI_SETUP: Final = """
mkdir -p $HOME/.icav2
echo "server-url: ica.illumina.com" > /root/.icav2/config.yaml

set +x
gcloud secrets versions access latest --secret=illumina_cpg_workbench_api --project=cpg-common | jq -r .apiKey > key
gcloud secrets versions access latest --secret=illumina_cpg_workbench_api --project=cpg-common | jq -r .projectID > projectID
echo "x-api-key: $(cat key)" >> $HOME/.icav2/config.yaml
icav2 projects enter $(cat projectID)
set -x
"""


coloredlogs.install(level=logging.INFO)


def calculate_needed_storage(
    cram: str,
) -> str:
    logging.info(f'Checking blob size for {cram}')
    storage_size: int = cpg_utils.to_path(cram).stat().st_size
    return f'{ceil(ceil((storage_size / (1024**3)) + 3) * 1.2)}Gi'


# No need to register this stage in Metamist I think, just ICA prep
@stage
class PrepareIcaForDragenAnalysis(SequencingGroupStage):
    """Set up ICA for a single realignment run.

    Creates a folder ID for the Dragen output to be written into.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> cpg_utils.Path:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return sg_bucket / GCP_FOLDER_FOR_ICA_PREP / f'{sequencing_group.name}_output_fid.json'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']

        prepare_ica_job: PythonJob = get_batch().new_python_job(
            name='PrepareIcaForDragenAnalysis',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        prepare_ica_job.image(image=image_path('ica'))

        outputs = self.expected_outputs(sequencing_group=sequencing_group)
        output_fids = prepare_ica_job.call(
            prepare_ica_for_analysis.run,
            ica_analysis_output_folder=config_retrieve(['ica', 'data_prep', 'output_folder']),
            api_root=ICA_REST_ENDPOINT,
            sg_name=sequencing_group.name,
            bucket_name=bucket_name,
        ).as_json()

        get_batch().write_output(
            output_fids,
            str(outputs),
        )
        return self.make_outputs(
            target=sequencing_group,
            data=outputs,
            jobs=prepare_ica_job,
        )


@stage
class UploadDataToIca(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> cpg_utils.Path:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return sg_bucket / GCP_FOLDER_FOR_ICA_PREP / f'{sequencing_group.name}_fids.json'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        upload_job = get_batch().new_bash_job(
            name='UploadDataToIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        upload_folder = config_retrieve(['ica', 'data_prep', 'upload_folder'])
        bucket: str = get_path_components_from_gcp_path(str(sequencing_group.cram))['bucket']

        upload_job.image(image=image_path('ica'))
        upload_job.storage(calculate_needed_storage(cram=str(sequencing_group.cram)))
        output = self.expected_outputs(sequencing_group=sequencing_group)
        authenticate_cloud_credentials_in_job(upload_job)

        # Check if the CRAM already exists in ICA before uploading. If it exists, just return the ID for the CRAM and CRAI
        upload_job.command(
            f"""
            {ICA_CLI_SETUP}
            cram_status=$(icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram --match-mode EXACT -o json | jq -r '.items[].id.details.status')
            if [[ $cram_status != "AVAILABLE" ]]
            then
                mkdir {sequencing_group.name}
                gcloud storage cp {sequencing_group.cram} .
                gcloud storage cp {sequencing_group.cram}.crai .
                icav2 projectdata upload {sequencing_group.name}.cram /{bucket}/{upload_folder}/{sequencing_group.name}/
                icav2 projectdata upload {sequencing_group.name}.cram.crai /{bucket}/{upload_folder}/{sequencing_group.name}/
            fi

            icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram --match-mode EXACT -o json | jq -r '.items[].id' > cram_id
            icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram.crai --match-mode EXACT -o json | jq -r '.items[].id' > crai_id

            jq -n --arg cram_id $(cat cram_id) --arg crai_id $(cat crai_id) '{{cram_fid: $cram_id, crai_fid: $crai_id}}' > {upload_job.ofile}
            """,
        )
        get_batch().write_output(upload_job.ofile, str(output))
        return self.make_outputs(
            target=sequencing_group,
            data=output,
            jobs=upload_job,
        )


@stage(required_stages=[PrepareIcaForDragenAnalysis, UploadDataToIca])
class AlignGenotypeWithDragen(SequencingGroupStage):
    # Output object with pipeline ID to GCP
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return sg_bucket / GCP_FOLDER_FOR_RUNNING_PIPELINE / f'{sequencing_group.name}_pipeline_id.json'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        dragen_pipeline_id = config_retrieve(['ica', 'pipelines', 'dragen_3_7_8'])
        dragen_ht_id: str = config_retrieve(['ica', 'pipelines', 'dragen_ht_id'])

        # Get the correct CRAM reference ID based off the choice made in the config
        cram_reference_id: str = config_retrieve(
            ['ica', 'cram_references', config_retrieve(['ica', 'cram_references', 'old_cram_reference'])],
        )

        user_tags: list[str] = config_retrieve(['ica', 'tags', 'user_tags'])
        technical_tags: list[str] = config_retrieve(['ica', 'tags', 'technical_tags'])
        reference_tags: list[str] = config_retrieve(['ica', 'tags', 'reference_tags'])
        user_reference: str = config_retrieve(['ica', 'tags', 'user_reference'])

        align_genotype_job: PythonJob = get_batch().new_python_job(
            name='AlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )
        align_genotype_job.image(image=image_path('ica'))
        outputs = self.expected_outputs(sequencing_group=sequencing_group)
        pipeline_call = align_genotype_job.call(
            run_align_genotype_with_dragen.run,
            ica_fids_path=inputs.as_path(target=sequencing_group, stage=UploadDataToIca),
            analysis_output_fid_path=inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis),
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            dragen_pipeline_id=dragen_pipeline_id,
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            api_root=ICA_REST_ENDPOINT,
        ).as_json()
        get_batch().write_output(
            pipeline_call,
            str(outputs),
        )
        return self.make_outputs(
            target=sequencing_group,
            data=outputs,
            jobs=align_genotype_job,
        )


@stage(analysis_type='dragen_align_genotype', required_stages=[PrepareIcaForDragenAnalysis, AlignGenotypeWithDragen])
class MonitorAlignGenotypeWithDragen(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> cpg_utils.Path:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return sg_bucket / GCP_FOLDER_FOR_RUNNING_PIPELINE / f'{sequencing_group.name}_pipeline_success.json'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        monitor_pipeline_run: PythonJob = get_batch().new_python_job(
            name='MonitorAlignGenotypeWithDragen',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )

        monitor_pipeline_run.image(image=image_path('ica'))
        outputs = self.expected_outputs(sequencing_group=sequencing_group)
        pipeline_run_results = monitor_pipeline_run.call(
            monitor_align_genotype_with_dragen.run,
            ica_pipeline_id_path=str(inputs.as_path(target=sequencing_group, stage=AlignGenotypeWithDragen)),
            api_root=ICA_REST_ENDPOINT,
        ).as_json()
        get_batch().write_output(
            pipeline_run_results,
            str(outputs),
        )
        return self.make_outputs(
            target=sequencing_group,
            data=outputs,
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


@stage(required_stages=[AlignGenotypeWithDragen])
class CancelIcaPipelineRun(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return (
            sg_bucket
            / GCP_FOLDER_FOR_RUNNING_PIPELINE
            / f'{sequencing_group.name}_pipeline_cancelled_at{slugify(str(datetime.now()))}.json'
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cancel_pipeline_run: PythonJob = get_batch().new_python_job(
            name='CancelIcaPipelineRun',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'Dragen'},
        )
        cancel_pipeline_run.image(image=image_path('ica'))
        outputs = self.expected_outputs(sequencing_group=sequencing_group)
        cancel_pipeline = cancel_pipeline_run.call(
            cancel_ica_pipeline_run.run,
            ica_pipeline_id_path=str(inputs.as_path(target=sequencing_group, stage=AlignGenotypeWithDragen)),
            api_root=ICA_REST_ENDPOINT,
        ).as_json()
        get_batch().write_output(
            cancel_pipeline,
            str(outputs),
        )
        return self.make_outputs(
            target=sequencing_group,
            data=outputs,
            jobs=cancel_pipeline_run,
        )


@stage(
    analysis_type='ica_data_download',
    required_stages=[PrepareIcaForDragenAnalysis, MonitorAlignGenotypeWithDragen, GvcfMlrWithDragen],
)
class DownloadDataFromIca(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        bucket_name: cpg_utils.Path = sequencing_group.dataset.prefix()
        return bucket_name / GCP_FOLDER_FOR_ICA_DOWNLOAD / f'{sequencing_group.name}'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']

        batch_instance: Batch = get_batch()
        ica_download_job: BashJob = batch_instance.new_bash_job(
            name='DownloadDataFromIca',
            attributes=(self.get_job_attrs() or {}) | {'tool': 'ICA'},
        )
        ica_analysis_folder_id_path: str = batch_instance.read_input(
            inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis),
        )
        ica_download_job.storage(storage=calculate_needed_storage(cram=str(sequencing_group.cram)))
        ica_download_job.memory('8Gi')
        ica_download_job.image(image=image_path('ica'))

        # Download an entire folder with ICA. Don't log projectId or API key
        authenticate_cloud_credentials_in_job(ica_download_job)
        ica_download_job.command(
            f"""
            {ICA_CLI_SETUP}
            icav2 projectdata download $(cat {ica_analysis_folder_id_path} | jq -r .analysis_output_fid) {sequencing_group.name} --exclude-source-path
            gcloud storage cp --recursive {sequencing_group.name} gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/{sequencing_group.name}
        """,
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=ica_download_job,
        )
