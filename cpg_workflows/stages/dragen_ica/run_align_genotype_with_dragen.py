import json
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

import cpg_utils
from cpg_workflows.stages.dragen_ica import ica_utils


def submit_dragen_run(
    cram_id: str,
    dragen_ht_id: str,
    cram_reference_id: str,
    qc_cross_cont_vcf_id: str,
    qc_cov_region_1_id: str,
    qc_cov_region_2_id: str,
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
                AnalysisDataInput(parameterCode='ref_tar', dataIds=[dragen_ht_id]),
                AnalysisDataInput(parameterCode='cram_reference', dataIds=[cram_reference_id]),
                AnalysisDataInput(parameterCode='qc_cross_cont_vcf', dataIds=[qc_cross_cont_vcf_id]),
                AnalysisDataInput(parameterCode='qc_coverage_region_1', dataIds=[qc_cov_region_1_id]),
                AnalysisDataInput(parameterCode='qc_coverage_region_2', dataIds=[qc_cov_region_2_id]),
            ],
            parameters=[
                AnalysisParameterInput(code='enable_map_align', value='true'),
                AnalysisParameterInput(code='enable_map_align_output', value='true'),
                AnalysisParameterInput(code='output_format', value='CRAM'),
                AnalysisParameterInput(code='enable_duplicate_marking', value='true'),
                AnalysisParameterInput(code='enable_variant_caller', value='true'),
                AnalysisParameterInput(code='vc_emit_ref_confidence', value='GVCF'),
                AnalysisParameterInput(code='vc_enable_vcf_output', value='false'),
                AnalysisParameterInput(code='enable_cnv', value='true'),
                AnalysisParameterInput(code='cnv_segmentation_mode', value='SLM'),
                AnalysisParameterInput(code='enable_sv', value='true'),
                AnalysisParameterInput(code='enable_cyp2d6', value='true'),
                AnalysisParameterInput(code='repeat_genotype_enable', value='false'),
                AnalysisParameterInput(
                    code='additional_args',
                    value="--read-trimmers polyg --soft-read-trimmers none --vc-hard-filter 'DRAGENHardQUAL:all:QUAL<5.0;LowDepth:all:DP<=1' --vc-frd-max-effective-depth 40 --vc-enable-joint-detection true --qc-coverage-ignore-overlaps true --qc-coverage-count-soft-clipped-bases true --qc-coverage-reports-1 cov_report,cov_report --qc-coverage-filters-1 'mapq<1,bq<0,mapq<1,bq<0'",
                ),
                AnalysisParameterInput(code='dragen_reports', value='false'),
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
    ica_fids_path: str,
    analysis_output_fid_path: str,
    dragen_ht_id: str,
    cram_reference_id: str,
    qc_cross_cont_vcf_id: str,
    qc_cov_region_1_id: str,
    qc_cov_region_2_id: str,
    dragen_pipeline_id: str,
    user_tags: list[str],
    technical_tags: list[str],
    reference_tags: list[str],
    user_reference: str,
    api_root: str,
    output_path: str,
) -> str:
    """_summary_

    Args:
        ica_fids_path (str): Path to the JSON in GCP holding the file IDs for the input CRAM and CRAI
        analysis_output_fid_path (str): Path to the JSON in GCP holding the folder ID to store the analysis outputs in ICA
        dragen_ht_id (str): The ICA file ID for the Dragen hash table used for mapping
        cram_reference_id (str): The ICA file ID for the FASTA reference that was used to align the CRAM file
        dragen_pipeline_id (str): The ICA pipeline ID for the Dragen pipeline that is going to be run
        user_tags (list[str]): List of user tags for the analysis (optional, can be empty)
        technical_tags (list[str]): List of technical tags for the analysis (optional, can be empty)
        reference_tags (list[str]): List of reference tags for the analysis (optional, can be empty)
        user_reference (str): A reference name for the pipeline run
        api_root (str): The ICA API root
        output_path (str): The path to write the pipeline ID to
    """

    SECRETS: dict[Literal['projectID', 'apiKey'], str] = ica_utils.get_ica_secrets()
    project_id: str = SECRETS['projectID']
    api_key: str = SECRETS['apiKey']

    coloredlogs.install(level=logging.INFO)
    configuration = icasdk.Configuration(host=api_root)
    configuration.api_key['ApiKeyAuth'] = api_key

    with open(cpg_utils.to_path(ica_fids_path), 'rt') as ica_fids_handle:
        ica_fids: dict[str, str] = json.load(ica_fids_handle)

    with open(cpg_utils.to_path(analysis_output_fid_path), 'rt') as analysis_outputs_fid_handle:
        analysis_output_fid: dict[str, str] = json.load(analysis_outputs_fid_handle)

    with icasdk.ApiClient(configuration=configuration) as api_client:
        api_instance = project_analysis_api.ProjectAnalysisApi(api_client)
        path_params: dict[str, str] = {'projectId': project_id}
        analysis_run_id: str = submit_dragen_run(
            cram_id=ica_fids['cram_fid'],
            dragen_ht_id=dragen_ht_id,
            cram_reference_id=cram_reference_id,
            qc_cross_cont_vcf_id=qc_cross_cont_vcf_id,
            qc_cov_region_1_id=qc_cov_region_1_id,
            qc_cov_region_2_id=qc_cov_region_2_id,
            dragen_pipeline_id=dragen_pipeline_id,
            ica_output_folder_id=analysis_output_fid['analysis_output_fid'],
            user_tags=user_tags,
            technical_tags=technical_tags,
            reference_tags=reference_tags,
            user_reference=user_reference,
            project_id=path_params,
            api_instance=api_instance,
        )

        logging.info(f'Submitted ICA run with pipeline ID: {analysis_run_id}')
        with cpg_utils.to_path(output_path).open('w') as f:
            f.write(analysis_run_id)

    return analysis_run_id
