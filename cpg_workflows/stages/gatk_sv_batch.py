"""
Driver for computing structural variants from GATK-SV from WGS data.
"""
from typing import Any

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path
from cpg_workflows.batch import get_batch
from cpg_workflows.stages.gatk_sv import (
    GatherSampleEvidence,
    make_combined_ped,
    get_ref_panel,
    get_images,
    get_references,
    add_gatk_sv_jobs,
)
from cpg_workflows.workflow import (
    stage,
    StageOutput,
    DatasetStage,
    StageInput,
    Dataset,
)

GATK_SV_COMMIT = 'b59a8b070da48ceed475814117787cf80cace170'
SV_CALLERS = ['manta', 'wham', 'scramble']
UNUSED_CALLERS = ['melt']


@stage(required_stages=[GatherSampleEvidence])
class GATKSVPipelineBatch(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        d: dict[str, Path] = {}

        fname_d = {
            'clean_vcf': f'clean.vcf.gz',
            'clean_vcf_index': f'clean.vcf.gz.tbi',
            'metrics_file_batch': f'metrics_file_batch.tsv',
            'qc_file': f'qc_file.tsv',
            'master_vcf_qc': f'master_vcf_qc.tar.gz',
            'final_sample_list': 'final_sample.list',
            'final_sample_outlier_list': 'final_sample_outlier.list',
        }

        for key, fname in fname_d.items():
            d[key] = self.prefix / dataset.name / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        sample_d = inputs.as_dict_by_target(GatherSampleEvidence)
        sids = dataset.get_sample_ids()

        input_dict: dict[str, Any] = {
            'name': dataset.name,
            'samples': dataset.get_sample_ids(),
            'counts_files_input': [
                str(sample_d[sid]['coverage_counts']) for sid in sids
            ],
            'sr_files_input': [str(sample_d[sid]['sr']) for sid in sids],
            'pe_files_input': [str(sample_d[sid]['pe']) for sid in sids],
            'sd_files_input': [str(sample_d[sid]['sd']) for sid in sids],
            'ped_file': str(make_combined_ped(dataset)),
            'chr_x': 'chrX',
            'chr_y': 'chrY',
            'MakeCohortVcf.min_sr_background_fail_batches': 0.5,
            'MakeCohortVcf.max_shard_size_resolve': 500,
            'MakeCohortVcf.max_shards_per_chrom_clean_vcf_step1': 200,
            'MakeCohortVcf.min_records_per_shard_clean_vcf_step1': 5000,
            'MakeCohortVcf.clean_vcf1b_records_per_shard': 10000,
            'MakeCohortVcf.samples_per_clean_vcf_step2_shard': 100,
            'MakeCohortVcf.clean_vcf5_records_per_shard': 5000,
            'MakeCohortVcf.random_seed': 0,
            'GATKSVPipelinePhase1.depth_interval_overlap': 0.8,
            'GATKSVPipelinePhase1.pesr_breakend_window': 300,
            'GATKSVPipelinePhase1.BAF_split_size': 10000,
            'GATKSVPipelinePhase1.RD_split_size': 10000,
            'GATKSVPipelinePhase1.PE_split_size': 10000,
            'GATKSVPipelinePhase1.SR_split_size': 1000,
            'GATKSVPipelinePhase1.depth_exclude_overlap_fraction': 0.5,
            'EvidenceQC.run_vcf_qc': True,
            'RegenotypeCNVs.n_RdTest_bins': 100000,
            'GATKSVPipelinePhase1.outlier_cutoff_nIQR': 999999,
            'GenotypeBatch.n_per_split': 5000,
            'GATKSVPipelinePhase1.ref_copy_number_autosomal_contigs': 2,
            'GenotypeBatch.n_RD_genotype_bins': 100000,
            'GATKSVPipelinePhase1.pesr_interval_overlap': 0.1,
            'GATKSVPipelinePhase1.gcnv_qs_cutoff': 30,
            'GATKSVPipelinePhase1.min_svsize': 50,
            'GATKSVPipelinePhase1.common_cnv_size_cutoff': 5000,
            'GATKSVPipelinePhase1.run_matrix_qc': True,
            'GATKSVPipelinePhase1.matrix_qc_distance': 1000000,
            'RegenotypeCNVs.n_per_split': 5000,
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs_input'] = [
                str(sample_d[sid][f'{caller}_vcf']) for sid in sids
            ]
            input_dict[f'use_{caller}'] = True
        for caller in UNUSED_CALLERS:
            input_dict[f'use_{caller}'] = False

        # reference panel gCNV models
        input_dict |= get_ref_panel(
            [
                'contig_ploidy_model_tar',
                'gcnv_model_tars',
                'qc_definitions',
            ]
        )
        input_dict |= get_images(
            [
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
                'cloud_sdk_docker',
            ]
        )
        input_dict |= get_references(
            [
                'genome_file',
                'primary_contigs_list',
                'primary_contigs_fai',
                'autosome_file',
                'allosome_file',
                'GATKSVPipelinePhase1.mei_bed',
                'GATKSVPipelinePhase1.rmsk',
                'GATKSVPipelinePhase1.cnmops_exclude_list',
                'GATKSVPipelinePhase1.segdups',
                'GATKSVPipelinePhase1.cytoband',
                {'GATKSVPipelinePhase1.pesr_exclude_intervals': 'pesr_exclude_list'},
                {'GATKSVPipelinePhase1.depth_exclude_intervals': 'depth_exclude_list'},
                'EvidenceQC.wgd_scoring_mask',
                'GatherSampleEvidenceBatch.wham_include_list_bed_file',
                'GatherSampleEvidenceBatch.preprocessed_intervals',
                'GatherSampleEvidenceBatch.manta_region_bed',
                {'GatherSampleEvidenceBatch.sd_locs_vcf': 'dbsnp_vcf'},
                'MakeCohortVcf.mei_bed',
                {'MakeCohortVcf.cytobands': 'cytoband'},
                'MakeCohortVcf.depth_exclude_list',
                {'MakeCohortVcf.pe_exclude_list': 'pesr_exclude_list'},
                'MakeCohortVcf.bin_exclude',
                'MakeCohortVcf.empty_file',
                'GenotypeBatch.bin_exclude',
            ]
        )
        ref_fasta = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
        input_dict |= {
            'reference_fasta': str(ref_fasta),
            'reference_index': str(ref_fasta) + '.fai',
            'reference_dict': str(ref_fasta.with_suffix('.dict')),
        }

        expected_d = self.expected_outputs(dataset)

        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)
