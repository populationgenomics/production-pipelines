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
from typing import Any

import click

from cpg_pipes import Path
from cpg_pipes.pipeline import stage, SampleStage, StageOutput, DatasetStage, \
    pipeline_click_options, create_pipeline
from cpg_pipes.pipeline.pipeline import StageInput, ExpectedResultT
from cpg_pipes.targets import Sample, Dataset, Cohort

from pipelines.gatk_sv.utils import add_gatksv_job, get_references, get_dockers, get_gcnv_models, SV_CALLERS

logger = logging.getLogger(__file__)


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
        """Add jobs to batch"""
        cram_path = sample.dataset.get_bucket() / 'cram' / f'{sample.id}.cram'

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(cram_path),
            'bam_or_cram_index': str(cram_path) + '.crai',
            'sample_id': sample.id,
            # 'revise_base_cram_to_bam': True,
        }

        input_dict |= get_dockers([
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
        ])
        input_dict |= get_references([
            'primary_contigs_fai',
            'primary_contigs_list',
            'reference_fasta',
            'reference_index',
            'reference_dict',
            'preprocessed_intervals',
            'manta_region_bed',
            'wham_include_list_bed_file',
        ])
        input_dict['reference_version'] = '38'

        expected_d = self.expected_outputs(sample)

        output_dict, j = add_gatksv_job(
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

        input_dict |= get_dockers([
            'sv_base_mini_docker',
            'sv_base_docker',
            'sv_pipeline_docker',
            'sv_pipeline_qc_docker',
        ])
        
        input_dict |= get_references([
            'genome_file',
            'wgd_scoring_mask',
        ])

        expected_d = self.expected_outputs(dataset)
        output_dict, j = add_gatksv_job(
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
            d[key] = dataset.get_bucket() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """Add jobs to batch"""
        d = inputs.as_dict_by_target(stage=GatherSampleEvidence)
        
        sids = dataset.get_sample_ids()
        
        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'ped_file': str(dataset.make_ped_file()),
            'counts': [str(d[sid]['coverage_counts']) for sid in sids],
            'PE_files': [str(d[sid]['pesr_disc']) for sid in sids],
            'SR_files': [str(d[sid]['pesr_split']) for sid in sids],
            # BAF generation
            # Required for cohorts if BAF_files not provided
            # Note: pipeline output is not sensitive to having some samples (~1%) missing BAF
            'gvcfs': [str(s.get_gvcf_path()) for s in dataset.get_samples()],
        }
        input_dict |= get_references([
            'genome_file',
            'primary_contigs_fai',
            'ref_fasta',
            'ref_fasta_index',
            'ref_dict',
            'inclusion_bed',
            'unpadded_intervals_file',
            'dbsnp_vcf',
            'dbsnp_vcf_index',
            'cnmops_chrom_file',
            'cnmops_exclude_list',
            'cnmops_allo_file',
            'cytoband',
            'mei_bed',
        ])
        input_dict |= get_gcnv_models()
        input_dict['ref_copy_number_autosomal_contigs'] = 2
        input_dict['allosomal_contigs'] = ['chrX', 'chrY'],
        input_dict['gcnv_qs_cutoff'] = 30,
        input_dict['min_svsize'] = 50,
        input_dict['matrix_qc_distance'] = 1000000,

        input_dict |= get_dockers([
            'sv_base_mini_docker',
            'sv_base_docker',
            'sv_pipeline_docker',
            'sv_pipeline_qc_docker',
            'linux_docker',
            'condense_counts_docker',
            'gatk_docker',
            'gcnv_gatk_docker',
            'cnmops_docker',
        ])

        expected_d = self.expected_outputs(dataset)
        output_dict, j = add_gatksv_job(
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
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """


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
