import logging
from typing import Any, Literal

import coloredlogs
import icasdk
from icasdk.apis.tags import project_analysis_api
from icasdk.model.analysis_data_input import AnalysisDataInput
from icasdk.model.analysis_parameter_input import AnalysisParameterInput
from icasdk.model.analysis_tag import AnalysisTag
from icasdk.model.create_nextflow_analysis import CreateNextflowAnalysis
from icasdk.model.nextflow_analysis_input import NextflowAnalysisInput

from cpg_workflows.stages.dragen_ica import ica_utils


def submit_dragen_run(
    cram_id: str,
    cram_index_id: str,
    dragen_ht_id: str,
    cram_reference_id: str,
    dragen_pipeline_id: str,
    ica_output_folder_id: str,
    user_tags: list[str],
    technical_tags: list[str],
    reference_tags: list[str],
    user_reference: str,
    project_id: dict[str, str],
    api_instance: project_analysis_api.ProjectAnalysisApi,
) -> str:
    header_params: dict[Any, Any] = {}
    body = CreateNextflowAnalysis(
        userReference=user_reference,
        pipelineId=dragen_pipeline_id,
        tags=AnalysisTag(
            technicalTags=technical_tags,
            userTags=user_tags,
            referenceTags=reference_tags,
        ),
        outputParentFolderId=ica_output_folder_id,
        analysisInput=NextflowAnalysisInput(
            inputs=[
                AnalysisDataInput(parameterCode='crams', dataIds=[cram_id]),
                AnalysisDataInput(parameterCode='crais', dataIds=[cram_index_id]),
                AnalysisDataInput(parameterCode='ref_tar', dataIds=[dragen_ht_id]),
                AnalysisDataInput(parameterCode='cram_reference', dataIds=[cram_reference_id]),
            ],
            parameters=[
                AnalysisParameterInput(
                    code='enable_map_align',
                    value='True',
                ),
                AnalysisParameterInput(code='enable_map_align_output', value='True'),
                AnalysisParameterInput(code='output_format', value='CRAM'),
                AnalysisParameterInput(code='enable_duplicate_marking', value='True'),
                AnalysisParameterInput(code='enable_variant_caller', value='True'),
                AnalysisParameterInput(code='vc_emit_reference_confidence', value='GVCF'),
                AnalysisParameterInput(code='vc_enable_vcf_output', value='False'),
                AnalysisParameterInput(code='enable_cnv', value='True'),
                AnalysisParameterInput(code='enable_sv', value='True'),
                AnalysisParameterInput(
                    code='additional_args',
                    value='--vc-hard-filter "DRAGENHardQUAL:all:QUAL<5.0;LowDepth:all:DP<=1" --vc-frd-max-effective-depth 40 --vc-enable-joint-detection true',
                ),
            ],
        ),
    )
    try:
        api_response = api_instance.create_nextflow_analysis(
            path_params=project_id,
            header_params=header_params,
            body=body,
        )
        return api_response.body['id']
    except icasdk.ApiException as e:
        raise icasdk.ApiException(f'Exception when calling ProjectAnalysisApi->create_nextflow_analysis: {e}') from e


def run(
    cram_name: str,
    cram_id_path: str,
    cram_index_id_path: str,
    dragen_ht_id: str,
    cram_reference_id: str,
    dragen_pipeline_id: str,
    ica_output_folder_id_path: str,
    user_tags: list[str],
    technical_tags: list[str],
    reference_tags: list[str],
    user_reference: str,
    gcp_bucket: str,
    pipeline_registration_path: str,
    api_root: str,
) -> None:
    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)
    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_analysis_api.ProjectAnalysisApi(api_client)
        path_params: dict[str, str] = {'projectId': project_id}
        analysis_run_id: str = submit_dragen_run(
            cram_id=ica_utils.read_blob_contents(full_blob_path=cram_id_path),
            cram_index_id=ica_utils.read_blob_contents(full_blob_path=cram_index_id_path),
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            dragen_pipeline_id=dragen_pipeline_id,
            ica_output_folder_id=ica_utils.read_blob_contents(full_blob_path=ica_output_folder_id_path),
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            project_id=path_params,
            api_instance=api_instance,
        )

        logging.info(f'Submitted ICA run with pipeline ID: {analysis_run_id}')
        ica_utils.register_output_to_gcp(
            bucket=gcp_bucket,
            object_contents=analysis_run_id,
            object_name=f'{cram_name}_pipeline_id',
            gcp_folder=pipeline_registration_path,
        )
