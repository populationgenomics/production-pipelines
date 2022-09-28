"""
Driver for computing structural variants from GATK-SV.
Uses metadata from the sample-metadata server.

Example running:

```
pipelines/gatk_sv.py 
    -n test \
    --analysis-dataset acute-care
    --dataset acute-care
```

OR with a config:

```
pipelines/gatk_sv.py --config pipelines/configs/acute-care-test.yml
```

WHERE pipelines/configs/acute-care-test.yml:

```
namespace: test
analysis_dataset: acute-care
datasets: [acute-care]
```
"""
import logging
import os
from os.path import join
from typing import Any

import click
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from cpg_utils import Path
from cpg_utils.config import set_config_paths, get_config
from cpg_utils.hail_batch import command, reference_path, image_path
from cpg_utils.workflows.batch import make_job_name, Batch
from cpg_utils.workflows.workflow import (
    stage,
    SampleStage,
    StageOutput,
    DatasetStage,
    CohortStage,
    StageInput,
    ExpectedResultT,
    Sample,
    Dataset,
    Cohort,
    get_workflow,
)

logger = logging.getLogger(__file__)


# GATK_SV_COMMIT = 'e4c5ffb1596e0de93eb491ae6d7cb6ec46874d11'
GATK_SV_COMMIT = 'fdb03a26c5f57f0619cd2d0bd2e9bc77550c175f'
SV_CALLERS = ['manta', 'wham']


def get_images(keys: list[str]) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.
    """
    return {k: image_path(k) for k in get_config()['images'].keys() if k in keys}


def get_gcnv_models() -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with gCNV models
    """
    res: dict[str, str | list[str]] = {}
    res['ref_panel_samples'] = get_config()['sv_ref_panel']['ref_panel_samples']
    res['contig_ploidy_model_tar'] = str(
        reference_path('broad/sv/ref_panel/contig_ploidy_model_tar')
    )
    res['gcnv_model_tars'] = [
        str(reference_path('broad/sv/ref_panel/model_tar_tmpl')).format(shard=i)
        for i in range(get_config()['sv_ref_panel']['model_tar_cnt'])
    ]
    res['ref_panel_PE_files'] = [
        str(reference_path('broad/sv/ref_panel/ref_panel_PE_file_tmpl')).format(
            sample=s
        )
        for s in get_config()['sv_ref_panel']['ref_panel_samples']
    ]
    res['ref_panel_SE_files'] = [
        str(reference_path('broad/sv/ref_panel/ref_panel_SE_file_tmpl')).format(
            sample=s
        )
        for s in get_config()['sv_ref_panel']['ref_panel_samples']
    ]
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
            if key == 'ref_fasta':
                ref_fasta = reference_path(f'broad/{ref_d_key}')
                res['reference_fasta'] = str(ref_fasta)
                res['reference_index'] = str(ref_fasta) + '.fai'
                res['reference_dict'] = str(ref_fasta.with_suffix('.dict'))
            else:
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
        Expected to produce coverage counts, and a VCF for each variant caller.
        """
        d: dict[str, Path] = dict()
        fname_by_key = {'coverage_counts': 'coverage_counts.tsv.gz'}
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
                sample.dataset.prefix()
                / 'gatk_sv'
                / stage_name
                / (sample.id + '-' + fname)
            )
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """Add jobs to batch"""
        cram_path = sample.dataset.prefix() / 'cram' / f'{sample.id}.cram'

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(cram_path),
            'bam_or_cram_index': str(cram_path) + '.crai',
            'sample_id': sample.id,
            # 'revise_base_cram_to_bam': True,
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
                'ref_fasta',
                'preprocessed_intervals',
                'manta_region_bed',
                'wham_include_list_bed_file',
                {'sd_locs_vcf': 'dbsnp_vcf'},
            ]
        )
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


@stage(required_stages=[GatherSampleEvidence, EvidenceQC])
class GatherBatchEvidence(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        d: dict[str, Path] = dict()
        fname_by_key = {
            'merged_BAF': '',
            'merged_BAF_index': '',
            'merged_SR': '',
            'merged_SR_index': '',
            'merged_PE': '',
            'merged_PE_index': '',
            'merged_bincov': '',
            'merged_bincov_index': '',
            'ploidy_matrix': '',
            'ploidy_plots': '',
            'combined_ped_file': '',
            'merged_dels': '',
            'merged_dups': '',
            'cnmops_del': '',
            'cnmops_del_index': '',
            'cnmops_dup': '',
            'cnmops_dup_index': '',
            'cnmops_large_del': '',
            'cnmops_large_del_index': '',
            'cnmops_large_dup': '',
            'cnmops_large_dup_index': '',
            'median_cov': '',
            'std_manta_vcf': '',
            'std_delly_vcf': '',
            'std_melt_vcf': '',
            'std_scramble_vcf': '',
            'std_wham_vcf': '',
            'PE_stats': '',
            'RD_stats': '',
            'SR_stats': '',
            'BAF_stats': '',
            'Matrix_QC_plot': '',
            'manta_tloc': '',
            'metrics_file_batchevidence': '',
        }
        for key, fname in fname_by_key.items():
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """Add jobs to batch"""
        d = inputs.as_dict_by_target(stage=GatherSampleEvidence)

        sids = dataset.get_sample_ids()

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'ped_file': str(
                dataset.write_ped_file(dataset.tmp_prefix() / 'samples.ped')
            ),
            'counts': [str(d[sid]['coverage_counts']) for sid in sids],
            'PE_files': [str(d[sid]['pesr_disc']) for sid in sids],
            'SR_files': [str(d[sid]['pesr_split']) for sid in sids],
        }
        input_dict |= get_references(
            [
                'genome_file',
                'primary_contigs_fai',
                'ref_fasta',
                'inclusion_bed',
                'unpadded_intervals_file',
                'dbsnp_vcf',
                'dbsnp_vcf_index',
                {'sd_locs_vcf': 'dbsnp_vcf'},
                {'cnmops_chrom_file': 'autosome_file'},
                'cnmops_exclude_list',
                {'cnmops_allo_file': 'allosome_file'},
                'cytobands',
                'mei_bed',
            ]
        )
        input_dict |= get_gcnv_models()
        input_dict['ref_copy_number_autosomal_contigs'] = 2
        input_dict['allosomal_contigs'] = ['chrX', 'chrY']
        input_dict['gcnv_qs_cutoff'] = 30
        input_dict['min_svsize'] = 50
        input_dict['matrix_qc_distance'] = 1000000

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
                'gcnv_gatk_docker',
                'cnmops_docker',
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


@stage(required_stages=GatherBatchEvidence)
class ClusterBatch(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#clusterbatch
    """

    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        Clustered SV VCFs
        Clustered depth-only call VCF
        """

        d: dict[str, Path] = dict()
        # TODO: Fill in fnames
        fname_by_key = {
            'clustered_depth_vcf' : '',
            'clustered_depth_vcf_index' : '',
            'clustered_manta_vcf' : '',
            'clustered_manta_vcf_index' : '',
            'clustered_wham_vcf' : '',
            'clustered_wham_vcf_index' : '',
            'clustered_melt_vcf' : '',
            'clustered_melt_vcf_index' : '',
            'clustered_scramble_vcf' : '',
            'clustered_scramble_vcf_index' : '',
            'metrics_file_clusterbatch' : '', 
        }
        for key, fname in fname_by_key.items():
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """

        # d = inputs.as_dict_by_target(GatherBatchEvidence)

        # sids = dataset.get_sample_ids()

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'ped_file': str(
                dataset.write_ped_file(dataset.tmp_prefix() / 'samples.ped')
            ),
        }

        # TODO: Define del_bed, del_dup  gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GATKSVPipelinePhase1/GATKSVPipelinePhase1/acce2c71-7458-4205-ae13-624f6efc9956/call-GatherBatchEvidence/GatherBatchEvidence/366e817a-feb2-4d3c-915c-a8bd25529b81/call-MergeDepth/MergeDepth/b58e04ee-846d-4be6-bba3-6dda2f2c995e/call-MergeSet_del/cacheCopy/ref_panel_1kg.DEL.bed.gz"
        # TODO: Define dup_bed "gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GATKSVPipelinePhase1/GATKSVPipelinePhase1/acce2c71-7458-4205-ae13-624f6efc9956/call-GatherBatchEvidence/GatherBatchEvidence/366e817a-feb2-4d3c-915c-a8bd25529b81/call-MergeDepth/MergeDepth/b58e04ee-846d-4be6-bba3-6dda2f2c995e/call-MergeSet_dup/cacheCopy/ref_panel_1kg.DUP.bed.gz",

        # TODO: Define wham_vcf_tar, manta_vcf_tar

        input_dict['depth_exclude_overlap_fraction'] = 0.5
        input_dict['depth_interval_overlap'] = 0.8
        input_dict['depth_clustering_algorithm'] = 'SINGLE_LINKAGE'
        input_dict['pesr_interval_overlap'] = 0.1
        input_dict['pesr_breakend_window'] = 300
        input_dict['pesr_clustering_algorithm'] = 'SINGLE_LINKAGE'

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_pipeline_docker',
                'linux_docker',
                'gatk_docker',
                'sv_pipeline_base_docker'
            ]
        )

        input_dict |= get_references(
            [
                'reference_fasta',  # ref_fasta
                'reference_fasta_fai',  # reference_index 
                'contig_list',  # primary_contigs_list
                'reference_dict',
                'pesr_exclude_intervals',  # pesr_exclude_list
                'depth_exclude_intervals',  # depth_exclude_list
                
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


@stage(required_stages=ClusterBatch)
class GenerateBatchMetrics(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#gatherbatchmetrics
    """
    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        Generates variant metrics for filtering.
        """
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Combined read count matrix, SR, PE, and BAF files (GatherBatchEvidence)
        Per-sample median coverage estimates (GatherBatchEvidence)
        Clustered SV VCFs (ClusterBatch)
        Clustered depth-only call VCF (ClusterBatch)
        """

        pass 


@stage(required_stages=GenerateBatchMetrics)
class FilterBatch(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#filterbatch
    """
    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        Random forest cutoffs file
        PED file with outlier samples excluded
        """
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Batch PED file
        Metrics file (GenerateBatchMetrics)
        Clustered SV and depth-only call VCFs (ClusterBatch)
        """

        pass 


@stage(required_stages=FilterBatch)
class MergeBatchSites(CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#mergebatchsites
    """
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        """
        Combined cohort PESR and depth VCFs
        """
        pass

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        List of filtered PESR VCFs (FilterBatch)
        List of filtered depth VCFs (FilterBatch)
        """

        pass 


@stage(required_stages=[GatherBatchEvidence, FilterBatch, MergeBatchSites])
class GenotypeBatch(CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#genotypebatch
    """
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        """
        Genotypes a batch of samples across unfiltered variants combined across all batches.
        Outputs: 
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        (Optional) Depth re-genotyping intervals list
        """
        pass

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Batch PESR and depth VCFs (FilterBatch)
        Cohort PESR and depth VCFs (MergeBatchSites)
        Batch read count, PE, and SR files (GatherBatchEvidence)
        """

        pass 


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
    get_workflow().run(stages=[GatherBatchEvidence])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
