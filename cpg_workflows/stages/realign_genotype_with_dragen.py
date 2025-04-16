import logging
from math import ceil
from typing import TYPE_CHECKING, Final

import coloredlogs

from hailtop.batch.job import BashJob

import cpg_utils
import cpg_utils.cloud
import cpg_utils.config
from cpg_utils.cloud import get_path_components_from_gcp_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import Batch, authenticate_cloud_credentials_in_job, command, get_batch
from cpg_workflows.stages.dragen_ica import (
    cancel_ica_pipeline_run,
    monitor_align_genotype_with_dragen,
    prepare_ica_for_analysis,
    run_align_genotype_with_dragen,
)
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.workflow import CohortStage, SequencingGroupStage, StageInput, StageOutput, stage

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


# @stage
# class CreateIcaPipelineFailureRecord(CohortStage):
#     # Note that this is only valid for cohorts from a single dataset.
#     def expected_outputs(self, cohort: Cohort) -> cpg_utils.Path:
#         sg_bucket: cpg_utils.Path = cohort.get_datasets()[0].prefix()
#         return sg_bucket / GCP_FOLDER_FOR_RUNNING_PIPELINE / f'All_ICA_failures_for_{cohort.name}.csv'

#     def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
#         batch_instance: Batch = get_batch()
#         ica_pipeline_failure_job: BashJob = batch_instance.new_bash_job(
#             name='CreateIcaPipelineFailureRecord',
#             attributes=(self.get_job_attrs(cohort) or {}) | {'tool': 'GCP'},
#         )

#         return self.make_outputs(
#             target=cohort,
#             data=self.expected_outputs(cohort=cohort),
#             jobs=ica_pipeline_failure_job,
#         )


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
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'ICA'},
        )
        prepare_ica_job.image(image=config_retrieve(['workflow', 'driver_image']))

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
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'ICA'},
        )
        upload_folder = config_retrieve(['ica', 'data_prep', 'upload_folder'])
        bucket: str = get_path_components_from_gcp_path(str(sequencing_group.cram))['bucket']

        upload_job.image(image=config_retrieve(['workflow', 'driver_image']))
        upload_job.storage(calculate_needed_storage(cram=str(sequencing_group.cram)))
        output = self.expected_outputs(sequencing_group=sequencing_group)
        authenticate_cloud_credentials_in_job(upload_job)

        # Check if the CRAM and CRAI already exists in ICA before uploading. If they exist, just return the ID for the CRAM and CRAI
        # The internal `command` method is a wrapper from cpg_utils.hail_batch that extends the normal hail batch command
        upload_job.command(
            command(
                f"""
                function copy_from_gcp {{
                    mkdir -p $BATCH_TMPDIR/{sequencing_group.name}
                    gcloud storage cp {sequencing_group.cram} $BATCH_TMPDIR/{sequencing_group.name}/
                    gcloud storage cp {sequencing_group.cram}.crai $BATCH_TMPDIR/{sequencing_group.name}/
                }}
                function upload_cram {{
                    icav2 projectdata upload $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram /{bucket}/{upload_folder}/{sequencing_group.name}/
                }}
                function upload_crai {{
                    icav2 projectdata upload $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.crai /{bucket}/{upload_folder}/{sequencing_group.name}/
                }}

                function get_fids {{
                    # Add a random delay before calling the ICA API to hopefully stop empty JSON files from being written to GCP
                    sleep $(shuf -i 1-30 -n 1)
                    icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram --match-mode EXACT -o json | jq -r '.items[].id' > cram_id
                    icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram.crai --match-mode EXACT -o json | jq -r '.items[].id' > crai_id

                    jq -n --arg cram_id $(cat cram_id) --arg crai_id $(cat crai_id) '{{cram_fid: $cram_id, crai_fid: $crai_id}}' > {upload_job.ofile}
                }}

                {ICA_CLI_SETUP}
                copy_from_gcp
                cram_status=$(icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram --match-mode EXACT -o json | jq -r '.items[].details.status')
                crai_status=$(icav2 projectdata list --parent-folder /{bucket}/{upload_folder}/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram.crai --match-mode EXACT -o json | jq -r '.items[].details.status')

                if [[ $cram_status != "AVAILABLE" ]]
                then
                    retry upload_cram
                fi

                if [[ $crai_status != "AVAILABLE" ]]
                then
                    retry upload_crai
                fi

                get_fids

                # Try 10 times to call the ICA API to get the required file ID data
                counter=0
                while [ $counter -le 10 ]
                    do
                    get_fids
                    if [ -s {upload_job.ofile} ]
                    then
                        break
                    elif [ $counter -eq 10 ]
                    then
                        exit 1
                    fi
                    counter=$((counter+1))
                done
                """,
                define_retry_function=True,
            ),
        )
        get_batch().write_output(upload_job.ofile, str(output))
        return self.make_outputs(
            target=sequencing_group,
            data=output,
            jobs=upload_job,
        )


@stage(
    required_stages=[PrepareIcaForDragenAnalysis, UploadDataToIca],
    analysis_type='dragen_align_genotype',
    analysis_keys=['pipeline_id'],
)
class ManageDragenPipeline(SequencingGroupStage):
    """
    Due to the nature of the Dragen pipeline and stage dependencies, we need to run, monitor and cancel the pipeline in the same stage.

    This stage handles the following tasks:
    1. Cancels a previous pipeline running on ICA if requested.
        - Set the `cancel_cohort_run` flag to `true` in the config and the stage will read the pipeline ID from the JSON file and cancel it.
    2. Resumes monitoring a previous pipeline run if it was interrupted.
        - Set the `monitor_previous` flag to `true` in the config. This will read the pipeline ID from the JSON file and monitor it.
    3. Initiates a new Dragen pipeline run if no previous run is found or if resuming is not requested.
    4. Monitors the progress of the Dragen pipeline run.
    """

    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> dict[str, cpg_utils.Path]:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()
        return {
            'output': sg_bucket / GCP_FOLDER_FOR_RUNNING_PIPELINE / f'{sequencing_group.name}_pipeline_success.json',
            'pipeline_id': sequencing_group.dataset.prefix()
            / GCP_FOLDER_FOR_RUNNING_PIPELINE
            / f'{sequencing_group.name}_pipeline_id.json',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        sg_bucket: cpg_utils.Path = sequencing_group.dataset.prefix()

        outputs = self.expected_outputs(sequencing_group=sequencing_group)

        stage_jobs: list[BashJob | PythonJob] = []

        if (
            config_retrieve(['ica', 'management', 'cancel_cohort_run'], False)
            and cpg_utils.to_path(outputs['pipeline_id']).exists()
        ):
            # Can only cancel a pipeline if the pipeline ID JSON exists
            with open(cpg_utils.to_path(outputs['pipeline_id']), 'rt') as pipeline_fid_handle:
                pipeline_id: str = pipeline_fid_handle.read().rstrip()
            # Cancel if requested
            logging.info('Cancelling pipeline run')
            cancel_job = get_batch().new_python_job(
                name='CancelIcaPipelineRun',
                attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'Dragen'},
            )
            cancel_job.image(image=config_retrieve(['workflow', 'driver_image']))
            cancel_pipeline_result = cancel_job.call(
                cancel_ica_pipeline_run.run,
                ica_pipeline_id_path=str(outputs['pipeline_id']),
                api_root=ICA_REST_ENDPOINT,
            ).as_json()
            get_batch().write_output(
                cancel_pipeline_result,
                str(
                    sg_bucket
                    / GCP_FOLDER_FOR_RUNNING_PIPELINE
                    / f'{sequencing_group.name}_pipeline_{pipeline_id}_cancelled.json',
                ),
            )
            stage_jobs.append(cancel_job)

            return self.make_outputs(
                target=sequencing_group,
                data=outputs,
                jobs=stage_jobs,
            )

        # Test if a previous pipeline should be re-monitored. Used for when monitor stage on batch crashes and we want to resume
        if (
            config_retrieve(['ica', 'management', 'monitor_previous'], False)
            and cpg_utils.to_path(
                outputs['pipeline_id'],
            ).exists()
        ):
            logging.info(f'Previous pipeline found for {sequencing_group.name}, not setting off a new one')
            with open(cpg_utils.to_path(outputs['pipeline_id']), 'rt') as pipeline_fid_handle:
                align_genotype_job_result: str = pipeline_fid_handle.read().rstrip()
        # We shouldn't actually want to catch this case, either we have a running pipeline in ICA with a pipeline ID registered in GCP, or we don't. We need both peices of information in order to resume.
        # elif resume:
        #     logging.warning(f'No previous pipeline found for {sequencing_group.name}, but resume flag set.')
        #     return self.make_outputs(target=sequencing_group)
        else:
            dragen_pipeline_id = config_retrieve(['ica', 'pipelines', 'dragen_3_7_8'])
            dragen_ht_id: str = config_retrieve(['ica', 'pipelines', 'dragen_ht_id'])

            # Get the correct CRAM reference ID based off the choice made in the config
            cram_reference_id: str = config_retrieve(
                ['ica', 'cram_references', config_retrieve(['ica', 'cram_references', 'old_cram_reference'])],
            )

            user_tags: list[str] = config_retrieve(['ica', 'tags', 'user_tags'])
            technical_tags: list[str] = config_retrieve(['ica', 'tags', 'technical_tags'])
            reference_tags: list[str] = config_retrieve(['ica', 'tags', 'reference_tags'])
            user_reference: str = sequencing_group.name

            qc_cross_cont_vcf_id = config_retrieve(['ica', 'qc', 'cross_cont_vcf'])
            qc_cov_region_1_id = config_retrieve(['ica', 'qc', 'coverage_region_1'])
            qc_cov_region_2_id = config_retrieve(['ica', 'qc', 'coverage_region_2'])

            align_genotype_job: PythonJob = get_batch().new_python_job(
                name='AlignGenotypeWithDragen',
                attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'Dragen'},
            )
            align_genotype_job.image(image=config_retrieve(['workflow', 'driver_image']))

            align_genotype_job_result = align_genotype_job.call(
                run_align_genotype_with_dragen.run,
                ica_fids_path=inputs.as_path(target=sequencing_group, stage=UploadDataToIca),
                analysis_output_fid_path=inputs.as_path(target=sequencing_group, stage=PrepareIcaForDragenAnalysis),
                dragen_ht_id=dragen_ht_id,
                cram_reference_id=cram_reference_id,
                qc_cross_cont_vcf_id=qc_cross_cont_vcf_id,
                qc_cov_region_1_id=qc_cov_region_1_id,
                qc_cov_region_2_id=qc_cov_region_2_id,
                dragen_pipeline_id=dragen_pipeline_id,
                user_tags=user_tags,
                technical_tags=technical_tags,
                reference_tags=reference_tags,
                user_reference=user_reference,
                api_root=ICA_REST_ENDPOINT,
                output_path=outputs['pipeline_id'],
            )
            stage_jobs.append(align_genotype_job)

        # now monitor that job
        monitor_pipeline_run: PythonJob = get_batch().new_python_job(
            name='MonitorAlignGenotypeWithDragen',
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'Dragen'},
        )

        monitor_pipeline_run.image(image=config_retrieve(['workflow', 'driver_image']))
        pipeline_run_results = monitor_pipeline_run.call(
            monitor_align_genotype_with_dragen.run,
            ica_pipeline_id=align_genotype_job_result,
            pipeline_id_file=str(outputs['pipeline_id']),
            api_root=ICA_REST_ENDPOINT,
        ).as_json()

        get_batch().write_output(
            pipeline_run_results,
            str(outputs['output']),
        )

        stage_jobs.append(monitor_pipeline_run)

        return self.make_outputs(
            target=sequencing_group,
            data=outputs,
            jobs=stage_jobs,
        )


@stage(analysis_type='dragen_mlr', required_stages=[ManageDragenPipeline])
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


@stage(
    analysis_type='cram',
    analysis_keys=['cram'],
    required_stages=[PrepareIcaForDragenAnalysis, ManageDragenPipeline],
)
class DownloadCramFromIca(SequencingGroupStage):
    """
    Download cram and crai files from ICA separately. This is to allow registrations of the cram files
    in metamist to be done via stage decorators. The pipeline ID needs to be read within the Hail BashJob to get the current
    pipeline ID. If read outside the job, it will get the pipeline ID from the previous pipeline run.
    """

    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        bucket_name: cpg_utils.Path = sequencing_group.dataset.prefix()
        return {
            'cram': bucket_name / GCP_FOLDER_FOR_ICA_DOWNLOAD / 'cram' / f'{sequencing_group.name}.cram',
            'crai': bucket_name / GCP_FOLDER_FOR_ICA_DOWNLOAD / 'cram' / f'{sequencing_group.name}.cram.crai',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']
        ica_analysis_output_folder = config_retrieve(['ica', 'data_prep', 'output_folder'])

        batch_instance: Batch = get_batch()
        ica_download_job: BashJob = batch_instance.new_bash_job(
            name='DownloadCramFromIca',
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'ICA'},
        )

        pipeline_id_path: cpg_utils.Path = inputs.as_path(
            target=sequencing_group,
            stage=ManageDragenPipeline,
            key='pipeline_id',
        )

        ica_download_job.storage(storage=calculate_needed_storage(cram=str(sequencing_group.cram)))
        ica_download_job.memory('8Gi')
        ica_download_job.image(image=config_retrieve(['workflow', 'driver_image']))

        # Download just the CRAM and CRAI files  with ICA. Don't log projectId or API key
        authenticate_cloud_credentials_in_job(ica_download_job)
        ica_download_job.command(
            command(
                f"""
                function download_cram {{
                cram_id=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram --match-mode EXACT -o json | jq -r '.items[].id')
                crai_id=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram.crai --match-mode EXACT -o json | jq -r '.items[].id')
                cram_md5=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.cram.md5sum --match-mode EXACT -o json | jq -r '.items[].id')
                icav2 projectdata download $cram_id $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram --exclude-source-path
                icav2 projectdata download $crai_id $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.crai --exclude-source-path
                icav2 projectdata download $cram_md5 $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.md5sum --exclude-source-path

                # Get md5sum of the downloaded CRAM file and compare it with the ICA md5sum
                # Checking here because using icav2 package to download which doesn't automatically perform checksum matching
                ica_md5_hash=$(cat $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.md5sum)
                cram_md5=$(cat $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram | md5sum | cut -d " " -f1)
                if [ "$cram_md5" != "$ica_md5_hash" ]; then
                    echo "Error: MD5 checksums do not match!"
                    echo "ICA MD5: $ica_md5_hash"
                    echo "Cram MD5: $cram_md5"
                    exit 1
                else
                    echo "MD5 checksums match."
                fi

                # Copy the CRAM and CRAI files to the bucket
                # Checksums are already checked by `gcloud storage cp`
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/cram/
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.crai gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/cram/
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.cram.md5sum gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/cram/
                }}

                {ICA_CLI_SETUP}
                mkdir -p $BATCH_TMPDIR/{sequencing_group.name}
                pipeline_id_filename=$(basename {pipeline_id_path})
                gcloud storage cp {pipeline_id_path} .
                pipeline_id=$(cat $pipeline_id_filename)
                echo "Pipeline ID: $pipeline_id"

                retry download_cram
                """,
                define_retry_function=True,
            ),
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=ica_download_job,
        )


@stage(
    analysis_type='gvcf',
    analysis_keys=['gvcf'],
    required_stages=[ManageDragenPipeline],
)
class DownloadGvcfFromIca(SequencingGroupStage):
    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        bucket_name: cpg_utils.Path = sequencing_group.dataset.prefix()
        return {
            'gvcf': bucket_name
            / GCP_FOLDER_FOR_ICA_DOWNLOAD
            / 'base_gvcf'
            / f'{sequencing_group.name}.hard-filtered.gvcf.gz',
            'gvcf_tbi': bucket_name
            / GCP_FOLDER_FOR_ICA_DOWNLOAD
            / 'base_gvcf'
            / f'{sequencing_group.name}.hard-filtered.gvcf.gz.tbi',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Download gVCF and gVCF TBI files from ICA separately. This is to allow registrations of the gVCF files
        in metamist to be done via stage decorators. The pipeline ID needs to be read within the Hail BashJob to get the current
        pipeline ID. If read outside the job, it will get the pipeline ID from the previous pipeline run.
        """
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']
        ica_analysis_output_folder = config_retrieve(['ica', 'data_prep', 'output_folder'])

        batch_instance: Batch = get_batch()
        ica_download_job: BashJob = batch_instance.new_bash_job(
            name='DownloadGvcfFromIca',
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'ICA'},
        )

        pipeline_id_path: cpg_utils.Path = inputs.as_path(
            target=sequencing_group,
            stage=ManageDragenPipeline,
            key='pipeline_id',
        )

        ica_download_job.storage(storage=calculate_needed_storage(cram=str(sequencing_group.cram)))
        ica_download_job.memory('8Gi')
        ica_download_job.image(image=config_retrieve(['workflow', 'driver_image']))

        # Download just the CRAM and CRAI files  with ICA. Don't log projectId or API key
        authenticate_cloud_credentials_in_job(ica_download_job)
        ica_download_job.command(
            command(
                f"""
                function download_gvcf {{
                gvcf_id=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.hard-filtered.gvcf.gz --match-mode EXACT -o json | jq -r '.items[].id')
                gvcf_tbi_id=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.hard-filtered.gvcf.gz.tbi --match-mode EXACT -o json | jq -r '.items[].id')
                gvcf_md5_id=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ --data-type FILE --file-name {sequencing_group.name}.hard-filtered.gvcf.gz.md5sum --match-mode EXACT -o json | jq -r '.items[].id')
                icav2 projectdata download $gvcf_id $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz --exclude-source-path
                icav2 projectdata download $gvcf_tbi_id $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz.tbi --exclude-source-path
                icav2 projectdata download $gvcf_md5_id $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz.md5sum --exclude-source-path

                # Get md5sum of the downloaded gVCF file and compare it with the ICA md5sum
                # Checking here because using icav2 package to download which doesn't automatically perform checksum matching
                ica_md5_hash=$(cat $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz.md5sum)
                gvcf_md5=$(cat $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz | md5sum | cut -d " " -f1)
                if [ "$gvcf_md5" != "$ica_md5_hash" ]; then
                    echo "Error: MD5 checksums do not match!"
                    echo "ICA MD5: $ica_md5_hash"
                    echo "Gvcf MD5: $gvcf_md5"
                    exit 1
                else
                    echo "MD5 checksums match."
                fi

                # Copy the gVCF and gVCF TBI files to the bucket
                # Checksums are already checked by `gcloud storage cp`
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/base_gvcf/
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz.tbi gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/base_gvcf/
                gcloud storage cp $BATCH_TMPDIR/{sequencing_group.name}/{sequencing_group.name}.hard-filtered.gvcf.gz.md5sum gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/base_gvcf/
                }}

                {ICA_CLI_SETUP}
                mkdir -p $BATCH_TMPDIR/{sequencing_group.name}
                pipeline_id_filename=$(basename {pipeline_id_path})
                gcloud storage cp {pipeline_id_path} .
                pipeline_id=$(cat $pipeline_id_filename)
                echo "Pipeline ID: $pipeline_id"

                retry download_gvcf
                """,
                define_retry_function=True,
            ),
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=ica_download_job,
        )


@stage(
    analysis_type='ica_data_download',
    required_stages=[ManageDragenPipeline, DownloadCramFromIca, DownloadGvcfFromIca],
)
class DownloadDataFromIca(SequencingGroupStage):
    """
    Download all files from ICA for a single realignment run except the CRAM and GVCF files.
    Register this batch download in Metamist.
    Does not register individual files in Metamist.
    """

    def expected_outputs(
        self,
        sequencing_group: SequencingGroup,
    ) -> cpg_utils.Path:
        bucket_name: cpg_utils.Path = sequencing_group.dataset.prefix()
        return bucket_name / GCP_FOLDER_FOR_ICA_DOWNLOAD / 'dragen_metrics' / f'{sequencing_group.name}'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        bucket_name: str = get_path_components_from_gcp_path(path=str(object=sequencing_group.cram))['bucket']
        ica_analysis_output_folder = config_retrieve(['ica', 'data_prep', 'output_folder'])
        pipeline_id_path: cpg_utils.Path = inputs.as_path(
            target=sequencing_group,
            stage=ManageDragenPipeline,
            key='pipeline_id',
        )

        batch_instance: Batch = get_batch()
        ica_download_job: BashJob = batch_instance.new_bash_job(
            name='DownloadDataFromIca',
            attributes=(self.get_job_attrs(sequencing_group) or {}) | {'tool': 'ICA'},
        )

        ica_download_job.storage(storage=calculate_needed_storage(cram=str(sequencing_group.cram)))
        ica_download_job.memory('8Gi')
        ica_download_job.image(image=config_retrieve(['workflow', 'driver_image']))

        # Download an entire folder (except crams and gvcfs) with ICA. Don't log projectId or API key
        authenticate_cloud_credentials_in_job(ica_download_job)
        ica_download_job.command(
            command(
                f"""
                function download_extra_data {{
                files_and_ids=$(icav2 projectdata list --parent-folder /{bucket_name}/{ica_analysis_output_folder}/{sequencing_group.name}/{sequencing_group.name}-$pipeline_id/{sequencing_group.name}/ -o json | jq -r '.items[] | select(.details.name | test(".cram|.gvcf") | not) | "\(.details.name) \(.id)"')
                while IFS= read -r line; do
                    name=$(echo "$line" | awk '{{print $1}}')
                    id=$(echo "$line" | awk '{{print $2}}')
                    echo "Downloading $name with ID $id"
                    icav2 projectdata download $id $BATCH_TMPDIR/{sequencing_group.name}/$name --exclude-source-path
                done <<< "$files_and_ids"
                gcloud storage cp --recursive $BATCH_TMPDIR/{sequencing_group.name}/* gs://{bucket_name}/{GCP_FOLDER_FOR_ICA_DOWNLOAD}/dragen_metrics/{sequencing_group.name}/
                }}

                {ICA_CLI_SETUP}
                # List all files in the folder except crams and gvcf and download them
                mkdir -p $BATCH_TMPDIR/{sequencing_group.name}
                pipeline_id_filename=$(basename {pipeline_id_path})
                gcloud storage cp {pipeline_id_path} .
                pipeline_id=$(cat $pipeline_id_filename)
                echo "Pipeline ID: $pipeline_id"

                retry download_extra_data
                """,
                define_retry_function=True,
            ),
        )

        return self.make_outputs(
            target=sequencing_group,
            data=self.expected_outputs(sequencing_group=sequencing_group),
            jobs=ica_download_job,
        )
