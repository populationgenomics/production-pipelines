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
from typing import Any

import click
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch import Batch

from cpg_pipes import images, Path
from cpg_pipes.hb.batch import make_job_name
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.pipeline import stage, SampleStage, StageOutput, DatasetStage, \
    CohortStage, pipeline_click_options, create_pipeline
from cpg_pipes.pipeline.pipeline import ExpectedResultT, StageInput
from cpg_pipes.targets import Sample, Dataset, Cohort

logger = logging.getLogger(__file__)


GATK_SV_COMMIT = 'c9e88f056fa154a01e2fcd7f8f0703342537a06e'
SV_CALLERS = ['manta', 'wham']

ACCESS_LEVEL = os.environ['CPG_ACCESS_LEVEL']


def get_dockers(keys: list[str]) -> dict[str, str]:
    """parse the WDL inputs with docker images"""
    with open(join(dirname(__file__), 'dockers.json')) as f:
        return {k: v for k, v in json.load(f).items() if k in keys}


def get_references(keys: list[str]) -> dict[str, str]:
    """parse the WDL inputs with reference files"""
    with open(join(dirname(__file__), 'references.json')) as f:
        return {k: v for k, v in json.load(f).items() if k in keys}


def add_gatksv_job(
    batch: Batch,
    dataset_name: str,
    wfl_name: str,
    input_dict: dict[str, Any], 
    expected_out_dict: dict[str, Path],
    sample_id: str | None = None,
):
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output. 
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset_name}'
    if sample_id: 
        output_prefix = join(output_prefix, sample_id)

    outputs_to_collect = dict()
    for key in expected_out_dict.keys():
        outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')
    
    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset_name)
    assert ACCESS_LEVEL
    output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=dataset_name,
        access_level=ACCESS_LEVEL,
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict={f'{wfl_name}.{k}': v for k, v in input_dict.items()},
        outputs_to_collect=outputs_to_collect,
        driver_image=images.DRIVER_IMAGE,
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
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
    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Expectes to produce coverage counts, and VCF for each variant caller
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'coverage_counts': 'coverage_counts.tsv.gz'
        }
        for caller in SV_CALLERS:
            fname_by_key[f'{caller}_vcf'] = f'{caller}.vcf.gz'
            fname_by_key[f'{caller}_index'] = f'{caller}.vcf.gz.tbi'

        # Split reads (SR) and Discordant read pairs (PE)
        for k in ['disc', 'split']:
            fname_by_key[f'pesr_{k}'] = f'pesr_{k}.txt.gz'
            fname_by_key[f'pesr_{k}_index'] = f'pesr_{k}.txt.gz.tbi'

        for key, fname in fname_by_key.items():
            stage_name = self.name.lower()
            d[key] = (
                sample.dataset.get_bucket() / 'gatk_sv' / stage_name /
                (sample.id + '-' + fname)
            )
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        cram_path = sample.dataset.get_bucket() / 'cram' / f'{sample.id}.cram'

        input_dict: dict[str, str] = {
            'bam_or_cram_file': str(cram_path),
            'bam_or_cram_index': str(cram_path) + '.crai',
            'sample_id': sample.id,
            # 'revise_base_cram_to_bam': True,
        }

        input_dict.update(get_dockers([
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
        ]))

        input_dict.update(get_references([
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
        ]))

        expected_d = self.expected_outputs(sample)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset_name=sample.dataset.name,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
        )

        return self.make_outputs(sample, data=output_dict, jobs=[j])


@stage(required_stages=GatherSampleEvidence)
class EvidenceQC(DatasetStage):
    """
    # https://github.com/broadinstitute/gatk-sv#evidenceqc
    """
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'ploidy_matrix':       'ploidy_matrix.bed.gz',
            'ploidy_plots':        'ploidy_plots.tar.gz',
            'WGD_dist':            'WGD_score_distributions.pdf',
            'WGD_matrix':          'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores':          'WGD_scores.txt.gz',
            'bincov_median':       'RD.txt.gz',
            'bincov_matrix':       'RD.txt.gz',
            'bincov_matrix_index': 'RD.txt.gz.tbi',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        for key, fname in fname_by_key.items():
            d[key] = dataset.get_bucket() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        d = inputs.as_dict_by_target(GatherSampleEvidence)

        sids = [s.id for s in dataset.get_samples()]

        input_dict = {
            'batch': dataset.name,
            'samples': sids,
            'run_vcf_qc': True,  # generates <caller>_qc_low/<caller>_qc_high
            'counts': [d[sid]['coverage_counts'] for sid in sids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [d[sid][f'{caller}_vcf'] for sid in sids]

        input_dict.update(get_dockers([
            'sv_base_mini_docker',
            'sv_base_docker',
            'sv_pipeline_docker',
            'sv_pipeline_qc_docker',
        ]))
        
        input_dict.update(get_references([
            'genome_file',
            'wgd_scoring_mask',
        ]))

        expected_d = self.expected_outputs(dataset)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset_name=dataset.name,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )

        return self.make_outputs(dataset, data=output_dict, jobs=[j])


@stage(required_stages=[GatherSampleEvidence])
class TrainGCNV(CohortStage):
    """
    # https://github.com/populationgenomics/gatk-sv/blob/main/wdl/TrainGCNV.wdl
    """
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to return 2 tarballs
        """
        d = dict()
        fname_by_key = {
            'cohort_contig_ploidy_model_tar': 'cohort_contig_ploidy_model.tar',
            'cohort_contig_ploidy_calls_tar': 'cohort_contig_ploidy_calls.tar',
        }
        for key, fname in fname_by_key.items():
            d[key] = cohort.analysis_dataset.get_bucket() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        d = inputs.as_dict_by_target(GatherSampleEvidence)
        
        sids = cohort.get_sample_ids()
        input_dict = {
            'cohort': cohort.target_id,
            'samples': sids,
            'count_files': [d[sid]['coverage_counts'] for sid in sids],
            'ref_copy_number_autosomal_contigs': 2,
            'num_intervals_per_scatter': 5000,
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [d[sid][f'{caller}_vcf'] for sid in sids]

        input_dict.update(get_dockers([
            'sv_base_mini_docker',
            'condense_counts_docker',
            'gatk_docker',
            'linux_docker',
        ]))
        
        input_dict.update(get_references([
            'reference_fasta',
            'reference_index',
            'reference_dict',
            'allosomal_contigs',
            'contig_ploidy_priors',
        ]))

        expected_d = self.expected_outputs(cohort)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset_name=cohort.analysis_dataset.name,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )

        return self.make_outputs(cohort, data=output_dict, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
): 
    """
    GATK-SV workflow.
    """
    pipeline = create_pipeline(
        name='gatk_sv',
        **kwargs
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
