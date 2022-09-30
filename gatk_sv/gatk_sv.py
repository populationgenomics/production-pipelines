"""
Driver for computing structural variants from GATK-SV from WGS data.
"""

import os
from os.path import join
from typing import Any

import click
import coloredlogs
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)

from cpg_utils import Path
from cpg_utils.config import set_config_paths, get_config
from cpg_utils.hail_batch import command, reference_path, image_path
from cpg_utils.workflows.batch import make_job_name, Batch
from cpg_utils.workflows.utils import exists
from cpg_utils.workflows.workflow import (
    stage,
    SampleStage,
    StageOutput,
    DatasetStage,
    StageInput,
    Sample,
    Dataset,
    run_workflow,
)

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='DEBUG', fmt=fmt)


GATK_SV_COMMIT = 'fdb03a26c5f57f0619cd2d0bd2e9bc77550c175f'
SV_CALLERS = ['manta', 'wham', 'scramble']


def get_images(keys: list[str]) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.
    """
    return {k: image_path(k) for k in get_config()['images'].keys() if k in keys}


def get_gcnv_models() -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with gCNV models
    """
    res: dict[str, str | list[str]] = {
        'ref_panel_samples': get_config()['sv_ref_panel']['ref_panel_samples'],
        'contig_ploidy_model_tar': str(
            reference_path('broad/sv/ref_panel/contig_ploidy_model_tar')
        ),
        'gcnv_model_tars': [
            str(reference_path('broad/sv/ref_panel/model_tar_tmpl')).format(shard=i)
            for i in range(get_config()['sv_ref_panel']['model_tar_cnt'])
        ],
        'ref_panel_PE_files': [
            str(reference_path('broad/sv/ref_panel/ref_panel_PE_file_tmpl')).format(
                sample=s
            )
            for s in get_config()['sv_ref_panel']['ref_panel_samples']
        ],
        'ref_panel_SE_files': [
            str(reference_path('broad/sv/ref_panel/ref_panel_SE_file_tmpl')).format(
                sample=s
            )
            for s in get_config()['sv_ref_panel']['ref_panel_samples']
        ],
    }
    return res


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        if isinstance(key, dict):
            key, ref_d_key = list(key.items())[0]
        else:
            ref_d_key = key
        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = str(reference_path(f'broad/sv/resources/{ref_d_key}'))
        except KeyError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))

    return res


def add_gatk_sv_job(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path],
    sample_id: str | None = None,
):
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sample_id:
        output_prefix = join(output_prefix, sample_id)

    outputs_to_collect = dict()
    for key in expected_out_dict.keys():
        outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset.name)
    output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=get_config()['workflow']['dataset'],
        access_level=get_config()['workflow']['access_level'],
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict={f'{wfl_name}.{k}': v for k, v in input_dict.items()},
        outputs_to_collect=outputs_to_collect,
        driver_image=image_path('hail'),
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(image_path('hail'))
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        cmds.append(f'gsutil cp $(cat {resource}) {out_path}')
    copy_j.command(command(cmds, setup_gcp=True))

    return output_dict, copy_j


@stage
class GatherSampleEvidence(SampleStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Expected to produce coverage counts, a VCF for each variant caller,
        and a txt for each type of SV evidence (SR, PE, SD).
        """
        d: dict[str, Path] = dict()

        # Coverage counts
        fname_by_key = {'coverage_counts': 'coverage_counts.tsv.gz'}

        # Caller's VCFs
        for caller in SV_CALLERS:
            fname_by_key[f'{caller}_vcf'] = f'{caller}.vcf.gz'
            fname_by_key[f'{caller}_index'] = f'{caller}.vcf.gz.tbi'

        # SV evidence
        # split reads:
        fname_by_key['pesr_split'] = 'sr.txt.gz'
        fname_by_key['pesr_split_index'] = 'sr.txt.gz.tbi'
        # discordant paired end reads:
        fname_by_key['pesr_disc'] = 'pe.txt.gz'
        fname_by_key['pesr_disc_index'] = 'pe.txt.gz.tbi'
        # site depth:
        fname_by_key['pesr_sd'] = 'sd.txt.gz'
        fname_by_key['pesr_sd_index'] = 'sd.txt.gz.tbi'

        for key, fname in fname_by_key.items():
            stage_name = self.name.lower()
            d[key] = (
                sample.dataset.prefix()
                / 'gatk_sv'
                / stage_name
                / (sample.id + '-' + fname)
            )
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """Add jobs to batch"""
        cram_path = sample.make_cram_path()
        assert exists(cram_path.path), cram_path

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(cram_path),
            'bam_or_cram_index': str(cram_path) + '.crai',
            'sample_id': sample.id,
        }

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'samtools_cloud_docker',
                'gatk_docker',
                'genomes_in_the_cloud_docker',
                'cloud_sdk_docker',
                'wham_docker',
                'manta_docker',
                'sv_pipeline_base_docker',
                'gatk_docker_pesr_override',
            ]
        )
        input_dict |= get_references(
            [
                'primary_contigs_fai',
                'primary_contigs_list',
                'preprocessed_intervals',
                'manta_region_bed',
                'wham_include_list_bed_file',
                {'sd_locs_vcf': 'dbsnp_vcf'},
            ]
        )
        input_dict |= {
            'reference_fasta': str(ref_fasta := reference_path(f'broad/ref_fasta')),
            'reference_index': str(ref_fasta) + '.fai',
            'reference_dict': str(ref_fasta.with_suffix('.dict')),
        }
        input_dict['reference_version'] = '38'

        expected_d = self.expected_outputs(sample)

        output_dict, j = add_gatk_sv_job(
            batch=self.b,
            dataset=sample.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
        )

        return self.make_outputs(sample, data=output_dict, jobs=[j])


@stage(required_stages=GatherSampleEvidence)
class EvidenceQC(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#evidenceqc
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'ploidy_matrix': 'ploidy_matrix.bed.gz',
            'ploidy_plots': 'ploidy_plots.tar.gz',
            'WGD_dist': 'WGD_score_distributions.pdf',
            'WGD_matrix': 'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores': 'WGD_scores.txt.gz',
            'bincov_median': 'RD.txt.gz',
            'bincov_matrix': 'RD.txt.gz',
            'bincov_matrix_index': 'RD.txt.gz.tbi',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        for key, fname in fname_by_key.items():
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """Add jobs to batch"""
        d = inputs.as_dict_by_target(GatherSampleEvidence)

        sids = dataset.get_sample_ids()

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'run_vcf_qc': True,  # generates <caller>_qc_low/<caller>_qc_high
            'counts': [str(d[sid]['coverage_counts']) for sid in sids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(d[sid][f'{caller}_vcf']) for sid in sids
            ]

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
            ]
        )

        input_dict |= get_references(
            [
                'genome_file',
                'wgd_scoring_mask',
            ]
        )

        expected_d = self.expected_outputs(dataset)
        output_dict, j = add_gatk_sv_job(
            batch=self.b,
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=output_dict, jobs=[j])


@click.command()
@click.argument('config_paths', nargs=-1)
def main(config_paths: list[str]):
    """
    Run a workflow, using CONFIG_PATHS in the order specified, overriding
    $CPG_CONFIG_PATH if specified.
    """
    if _cpg_config_path_env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths = _cpg_config_path_env_var.split(',') + config_paths
    set_config_paths(list(config_paths))
    run_workflow(stages=[EvidenceQC])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
