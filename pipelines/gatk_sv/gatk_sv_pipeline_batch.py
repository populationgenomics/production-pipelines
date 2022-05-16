"""
Driver for computing structural variants from GATK-SV.
Uses metadata from the sample-metadata server.
Runs the entire pipeline in Batch mode with GATKSVPipelineBatch.py.

Example running:

```
pipelines/gatk_sv_pipeline_batch.py
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
from typing import Any

import click

from cpg_pipes import Path
from cpg_pipes.pipeline import stage, StageOutput, DatasetStage, \
    pipeline_click_options, create_pipeline
from cpg_pipes.pipeline.pipeline import StageInput
from cpg_pipes.targets import Dataset
from pipelines.gatk_sv.utils import add_gatksv_job, get_references, get_dockers, \
    get_gcnv_models, SV_CALLERS

logger = logging.getLogger(__file__)


@stage
class GATKSVPipelineBatch(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GATKSVPipelineBatch.wdl
    """
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        File clean_vcf = MakeCohortVcf.vcf
        File clean_vcf_index = MakeCohortVcf.vcf_index
        File metrics_file_batch = CatBatchMetrics.out
        File qc_file = BatchQC.out
        File master_vcf_qc = MakeCohortVcf.vcf_qc
        File? metrics_file_makecohortvcf = MakeCohortVcf.metrics_file_makecohortvcf
        File final_sample_list = GATKSVPipelinePhase1.batch_samples_postOutlierExclusion_file
        File final_sample_outlier_list = GATKSVPipelinePhase1.outlier_samples_excluded_file
    
        # Additional outputs for creating a reference panel
        Array[File] counts = counts_files_
        Array[File] PE_files = pe_files_
        Array[File] PE_files_index = pe_files_index_
        Array[File] SR_files = sr_files_
        Array[File] SR_files_index = sr_files_index_
        Array[File]? manta_vcfs = manta_vcfs_
        Array[File]? manta_vcfs_index = manta_vcfs_index_
        Array[File]? melt_vcfs = melt_vcfs_
        Array[File]? melt_vcfs_index = melt_vcfs_index_
        Array[File]? wham_vcfs = wham_vcfs_
        Array[File]? wham_vcfs_index = wham_vcfs_index_
    
        File medianfile = GATKSVPipelinePhase1.median_cov
        File merged_coverage_file = GATKSVPipelinePhase1.merged_bincov
        File merged_coverage_file_index = GATKSVPipelinePhase1.merged_bincov_index
        File merged_baf_file = GATKSVPipelinePhase1.merged_BAF
        File merged_baf_file_index = GATKSVPipelinePhase1.merged_BAF_index
        File merged_disc_file = GATKSVPipelinePhase1.merged_PE
        File merged_disc_file_index = GATKSVPipelinePhase1.merged_PE_index
        File merged_split_file = GATKSVPipelinePhase1.merged_SR
        File merged_split_file_index = GATKSVPipelinePhase1.merged_SR_index
    
        File del_bed = GATKSVPipelinePhase1.merged_dels
        File del_bed_index = GATKSVPipelinePhase1.merged_dels + ".tbi"
        File dup_bed = GATKSVPipelinePhase1.merged_dups
        File dup_bed_index = GATKSVPipelinePhase1.merged_dups + ".tbi"
        Array[File]? std_manta_vcfs = GATKSVPipelinePhase1.std_manta_vcf
        Array[File]? std_melt_vcfs = GATKSVPipelinePhase1.std_melt_vcf
        Array[File]? std_scramble_vcfs = GATKSVPipelinePhase1.std_scramble_vcf
        Array[File]? std_wham_vcfs = GATKSVPipelinePhase1.std_wham_vcf
    
        File merged_depth_vcf = GATKSVPipelinePhase1.depth_vcf
        File merged_depth_vcf_index = GATKSVPipelinePhase1.depth_vcf_index
        File? merged_manta_vcf = GATKSVPipelinePhase1.manta_vcf
        File? merged_manta_vcf_index = GATKSVPipelinePhase1.manta_vcf_index
        File? merged_melt_vcf = GATKSVPipelinePhase1.melt_vcf
        File? merged_melt_vcf_index = GATKSVPipelinePhase1.melt_vcf_index
        File? merged_wham_vcf = GATKSVPipelinePhase1.wham_vcf
        File? merged_wham_vcf_index = GATKSVPipelinePhase1.wham_vcf_index
    
        File evidence_metrics = GATKSVPipelinePhase1.evidence_metrics
        File evidence_metrics_common = GATKSVPipelinePhase1.evidence_metrics_common
    
        File filtered_depth_vcf = select_first([GATKSVPipelinePhase1.filtered_depth_vcf])
        File filtered_pesr_vcf = select_first([GATKSVPipelinePhase1.filtered_pesr_vcf])
        File cohort_pesr_vcf = select_first([GATKSVPipelinePhase1.filtered_pesr_vcf])
        File cohort_depth_vcf = select_first([GATKSVPipelinePhase1.filtered_depth_vcf])
        File? sites_filtered_manta_vcf = GATKSVPipelinePhase1.sites_filtered_manta_vcf
        File? sites_filtered_wham_vcf = GATKSVPipelinePhase1.sites_filtered_wham_vcf
        File? sites_filtered_melt_vcf = GATKSVPipelinePhase1.sites_filtered_melt_vcf
        File? sites_filtered_depth_vcf = GATKSVPipelinePhase1.sites_filtered_depth_vcf
    
        File cutoffs = GATKSVPipelinePhase1.cutoffs
        File genotyped_pesr_vcf = GenotypeBatch.genotyped_pesr_vcf
        File genotyped_depth_vcf = GenotypeBatch.genotyped_depth_vcf
        File regeno_coverage_medians = GenotypeBatch.regeno_coverage_medians
        File regenotyped_depth_vcf = RegenotypeCNVs.regenotyped_depth_vcfs[0]
    
        File genotype_pesr_pesr_sepcutoff = select_first([GenotypeBatch.trained_genotype_pesr_pesr_sepcutoff])
        File genotype_pesr_depth_sepcutoff = select_first([GenotypeBatch.trained_genotype_pesr_depth_sepcutoff])
        File genotype_depth_pesr_sepcutoff = select_first([GenotypeBatch.trained_genotype_depth_pesr_sepcutoff])
        File genotype_depth_depth_sepcutoff = select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])
        File depth_gt_rd_sep_file = select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])
        File PE_metrics = select_first([GenotypeBatch.trained_PE_metrics])
        File SR_metrics = select_first([GenotypeBatch.trained_SR_metrics])
        File raw_sr_bothside_pass_file = GenotypeBatch.sr_bothside_pass
        File raw_sr_background_fail_file = GenotypeBatch.sr_background_fail        
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'clean_vcf': '',
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
        """Add jobs to batch"""
        sids = dataset.get_sample_ids()
        crams = [s.get_cram_path() for s in dataset.get_samples()]

        input_dict: dict[str, Any] = {}

        input_dict |= get_references([
            'genome_file',
            'primary_contigs_list',
            'primary_contigs_fai',
            'reference_fasta',
            'reference_index',
            'reference_dict',
            'autosome_file',
            'allosome_file',
            'qc_definitions',
            'GATKSVPipelinePhase1.unpadded_intervals_file',
            'GATKSVPipelinePhase1.dbsnp_vcf',
            'GATKSVPipelinePhase1.dbsnp_vcf_index',
            'GATKSVPipelinePhase1.cnmops_exclude_list',
            'GATKSVPipelinePhase1.rmsk',
            'GATKSVPipelinePhase1.pesr_exclude_list',
            'GATKSVPipelinePhase1.mei_bed',
            'GATKSVPipelinePhase1.inclusion_bed',
            'GATKSVPipelinePhase1.cytoband',
            'GATKSVPipelinePhase1.segdups',
            'GatherSampleEvidenceBatch.preprocessed_intervals',
            'GatherSampleEvidenceBatch.wham_include_list_bed_file',
            'GatherSampleEvidenceBatch.manta_region_bed',
            'EvidenceQC.wgd_scoring_mask',
            'GenotypeBatch.bin_exclude',
            'MakeCohortVcf.depth_exclude_list',
            {'MakeCohortVcf.pe_exclude_list': 'pesr_exclude_list'},
            'MakeCohortVcf.bin_exclude',
            'MakeCohortVcf.mei_bed',
            {'MakeCohortVcf.cytobands': 'cytoband'},
            'MakeCohortVcf.empty_file',
        ])

        gcnv_models_d = get_gcnv_models()
        input_dict |= {
            'contig_ploidy_model_tar': str(gcnv_models_d['contig_ploidy_model_tar']),
            'gcnv_model_tars': gcnv_models_d['gcnv_model_tars'],
        }
        
        input_dict |= {
            'name': dataset.name,
            'samples': sids,
            'bam_or_cram_files': [str(cram.path) for cram in crams],
            'bam_or_cram_indexes': [str(cram.index_path) for cram in crams],
            'use_manta': True,
            'use_melt': False,
            'use_scramble': False,
            'use_wham': True,
            'gvcfs': [str(s.get_gvcf_path().path) for s in dataset.get_samples()],
            'ped_file': str(dataset.make_ped_file()),
        }
        # All values below are for some reason not specified in the high-level 
        # workflow "input", yet are required, and I've no idea what any of those 
        # mean. So I just pasted values from GATK-SV tests and examples I found
        # in the repo.
        input_dict |= {
            'MakeCohortVcf.clean_vcf1b_records_per_shard': 10000,
            'MakeCohortVcf.chr_y': 'chrY',
            'MakeCohortVcf.chr_x': 'chrX',
            'MakeCohortVcf.max_shards_per_chrom_clean_vcf_step1': 200,
            'MakeCohortVcf.max_shard_size_resolve': 500,
            'MakeCohortVcf.clean_vcf5_records_per_shard': 5000,
            'MakeCohortVcf.min_sr_background_fail_batches': 0.5,
            'MakeCohortVcf.samples_per_clean_vcf_step2_shard': 100,
            'MakeCohortVcf.min_records_per_shard_clean_vcf_step1': 5000,
            'GenotypeBatch.n_RD_genotype_bins': 100000,
            'GenotypeBatch.n_per_split': 500,
            'GATKSVPipelinePhase1.common_cnv_size_cutoff': 5000,
            'GATKSVPipelinePhase1.run_matrix_qc': False,
            'GATKSVPipelinePhase1.depth_frac': 0.8,
            'GATKSVPipelinePhase1.pesr_frac': 0.1,
            'GATKSVPipelinePhase1.BAF_split_size': 10000,
            'GATKSVPipelinePhase1.depth_flags': '--merge-coordinates',
            'GATKSVPipelinePhase1.min_svsize': 50,
            'GATKSVPipelinePhase1.pesr_flags': '--preserve-ids',
            'GATKSVPipelinePhase1.outlier_cutoff_nIQR': 10000,
            'GATKSVPipelinePhase1.gcnv_qs_cutoff': 30,
            'GATKSVPipelinePhase1.pesr_distance': 300,
            'GATKSVPipelinePhase1.PE_split_size': 10000,
            'GATKSVPipelinePhase1.SR_split_size': 1000,
            'GATKSVPipelinePhase1.RD_split_size': 10000,
            'GATKSVPipelinePhase1.pesr_svsize': 0,
            'GATKSVPipelinePhase1.ref_copy_number_autosomal_contigs': 2,
            'GATKSVPipelinePhase1.matrix_qc_distance': 1000000,
            'EvidenceQC.run_vcf_qc': False,
            'RegenotypeCNVs.n_RdTest_bins': 1,
            'RegenotypeCNVs.n_per_split': 5000,
        }

        input_dict |= get_dockers([
            'sv_base_mini_docker',
            'sv_base_docker',
            'sv_pipeline_docker',
            'sv_pipeline_hail_docker',
            'sv_pipeline_updates_docker',
            'sv_pipeline_rdtest_docker',
            'sv_pipeline_base_docker',
            'sv_pipeline_qc_docker',
            'linux_docker',
            'cnmops_docker',
            'gatk_docker',
            'condense_counts_docker',
            'genomes_in_the_cloud_docker',
            'samtools_cloud_docker',
            'manta_docker',
            'scramble_docker',
            'wham_docker',
            'cloud_sdk_docker',
        ])

        expected_d: dict[str, Path] = self.expected_outputs(dataset)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=output_dict, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
): 
    """
    GATK-SV workflow.
    """
    kwargs['version'] = 'v0-1'
    pipeline = create_pipeline(
        name='gatk_sv',
        **kwargs
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
