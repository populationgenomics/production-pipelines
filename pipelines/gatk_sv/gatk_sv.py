"""
Driver for computing structural variants from GATK-SV.
Uses metadata from the sample-metadata server.

Example running:
gatk_sv.py \
    -n test \
    --analysis-project seqr \
    --input-project acute-care \
    -s CPG199869 -s CPG199877 -s CPG199885 \
    -s CPG199893 -s CPG199901 -s CPG199919 \
    -s CPG199927 -s CPG199935 -s CPG199943 \
    -s CPG54098
"""

import json
import logging
import os
from os.path import join, dirname
from typing import Dict, List, Optional

import click
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)

from cpg_pipes import images
from cpg_pipes.hb.batch import job_name
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.pipeline import (
    Sample,
    SampleStage,
    pipeline_click_options,
    StageInput,
    stage,
    Pipeline,
    find_stages_in_module,
    StageOutput,
)

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


GATK_SV_COMMIT = 'c9e88f056fa154a01e2fcd7f8f0703342537a06e'
SV_CALLERS = ['manta', 'wham']

ACCESS_LEVEL = os.environ['ACCESS_LEVEL']


def get_dockers(keys: List[str]) -> Dict[str, object]:
    """parse the WDL inputs with docker images"""
    with open(join(dirname(__file__), 'dockers.json')) as f:
        return {k: v for k, v in json.load(f).items() if k in keys}


def get_references(keys: List[str]) -> Dict[str, object]:
    """parse the WDL inputs with reference files"""
    with open(join(dirname(__file__), 'references.json')) as f:
        return {k: v for k, v in json.load(f).items() if k in keys}


def add_gatksv_job(
    pipeline,
    wfl_name: str,
    input_dict: Dict[str, object],
    expected_out_dict: Dict[str, str],
    project_name: Optional[str] = None,
    sample_id: Optional[str] = None,
):
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_suffix = f'gatk_sv/output/{wfl_name}'
    if project_name:
        output_suffix = join(output_suffix, project_name)
    if sample_id:
        output_suffix = join(output_suffix, sample_id)

    outputs_to_collect = dict()
    for key in expected_out_dict.keys():
        outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

    job_prefix = job_name(wfl_name, sample=sample_id, project=project_name)
    output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=pipeline.b,
        job_prefix=job_prefix,
        dataset=pipeline.analysis_project.stack,
        access_level=ACCESS_LEVEL,
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_suffix=output_suffix,
        input_dict={f'{wfl_name}.{k}': v for k, v in input_dict.items()},
        outputs_to_collect=outputs_to_collect,
        driver_image=images.DRIVER_IMAGE,
    )

    copy_j = pipeline.b.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(images.DRIVER_IMAGE)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        cmds.append(f'gsutil cp $(cat {resource}) {out_path}')
    copy_j.command(wrap_command(cmds, setup_gcp=True))

    return output_dict, copy_j


@stage
class GatherSampleEvidence(SampleStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    """

    def expected_result(self, sample: Sample) -> Dict[str, str]:
        d = dict()
        fname_by_key = {'coverage_counts': 'coverage_counts.tsv.gz'}
        for caller in SV_CALLERS:
            fname_by_key[f'{caller}_vcf'] = f'{caller}.vcf.gz'
            fname_by_key[f'{caller}_index'] = f'{caller}.vcf.gz.tbi'

        # Split reads (SR) and Discordant read pairs (PE)
        for k in ['disc', 'split']:
            fname_by_key[f'pesr_{k}'] = f'pesr_{k}.txt.gz'
            fname_by_key[f'pesr_{k}_index'] = f'pesr_{k}.txt.gz.tbi'

        for key, fname in fname_by_key.items():
            stage = self.name.lower()
            d[key] = join(
                sample.project.get_bucket(), 'gatk_sv', stage, sample.id + '-' + fname
            )
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = f'{sample.project.get_bucket()}/cram/{sample.id}.cram'

        input_dict: Dict[str, object] = {
            'bam_or_cram_file': cram_path,
            'bam_or_cram_index': cram_path + '.crai',
            'sample_id': sample.id,
            # 'revise_base_cram_to_bam': True,
        }

        input_dict.update(
            get_dockers(
                [
                    'sv_pipeline_docker',
                    'sv_base_mini_docker',
                    'samtools_cloud_docker',
                    'gatk_docker',
                    'genomes_in_the_cloud_docker',
                    'cloud_sdk_docker',
                    'wham_docker',
                    # 'melt_docker',
                    'manta_docker',
                    'sv_pipeline_base_docker',
                    'gatk_docker_pesr_override',
                ]
            )
        )

        input_dict.update(
            get_references(
                [
                    'primary_contigs_fai',
                    'primary_contigs_list',
                    'reference_fasta',
                    'reference_index',
                    'reference_dict',
                    'reference_version',
                    'preprocessed_intervals',
                    'manta_region_bed',
                    'wham_include_list_bed_file',
                    # 'melt_standard_vcf_header',
                ]
            )
        )

        expected_d = self.expected_result(sample)

        output_dict, j = add_gatksv_job(
            self.pipe,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
            project_name=sample.project.name,
        )

        return self.make_outputs(sample, data=output_dict, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    input_projects: List[str],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_projects
    title = (
        f'GATK-SV: samples from: {", ".join(input_projects)}'
        f', version {output_version}'
    )
    pipeline = Pipeline(
        name='gatk_sv',
        title=title,
        input_projects=input_projects,
        output_version=output_version,
        **kwargs,
    )
    pipeline.set_stages(find_stages_in_module(__name__))
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
