"""
The first multi-sample workflow, containing stages from
GenerateBatchMetrics to Genotype Batch

Though, there is a separate process required to merge files
across batches between FilterBatch and GenotypeBatch
"""


from typing import Any

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import stage, StageOutput, StageInput, Cohort, CohortStage

from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
    get_ref_panel,
    make_combined_ped,
    SV_CALLERS,
)


@stage
class GatherBatchEvidence(CohortStage):
    """
    This requires restriction by Samples, and runs separately from
    the initial variant calling and EvidenceQC

    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """create the output paths for GatherBatchEvidence"""
        ending_by_key = {
            'cnmops_dup': 'DUP.header.bed.gz',
            'cnmops_dup_index': 'DUP.header.bed.gz.tbi',
            'cnmops_del': 'DEL.header.bed.gz',
            'cnmops_del_index': 'DEL.header.bed.gz.tbi',
            'cnmops_large_del': 'DEL.large.bed.gz',
            'cnmops_large_del_index': 'DEL.large.bed.gz.tb',
            'cnmops_large_dup': 'DUP.large.bed.gz',
            'cnmops_large_dup_index': 'DUP.large.bed.gz.tbi',
            'merged_SR': 'sr.txt.gz',
            'merged_SR_index': 'sr.txt.gz.tbi',
            'merged_PE': 'pe.txt.gz',
            'merged_PE_index': 'pe.txt.gz.tbi',
            'merged_BAF': 'baf.txt.gz',
            'merged_BAF_index': 'baf.txt.gz.tbi',
            'merged_bincov': 'RD.txt.gz',
            'merged_bincov_index': 'RD.txt.gz.tbi',
            'SR_stats': 'SR.QC_matrix.txt',
            'PE_stats': 'PE.QC_matrix.txt',
            'BAF_stats': 'BAF.QC_matrix.txt',
            'RD_stats': 'RD.QC_matrix.txt',
            'median_cov': 'medianCov.transposed.bed',
            'merged_dels': 'DEL.bed.gz',
            'merged_dups': 'DUP.bed.gz',
            'Matrix_QC_plot': '00_matrix_FC_QC.png',
        }
        for caller in SV_CALLERS:
            ending_by_key[f'std_{caller}_vcf_tar'] = f'{caller}.tar.gz'

        return {key: self.prefix / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """Add jobs to Batch"""
        samples = cohort.get_samples(only_active=True)

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'samples': [sam.id for sam in samples],
            'ped_file': str(make_combined_ped(cohort, self.prefix)),
            'counts': [
                str(
                    sample.make_sv_evidence_path / f'{sample.id}.coverage_counts.tsv.gz'
                )
                for sample in samples
            ],
            'SR_files': [
                str(sample.make_sv_evidence_path / f'{sample.id}.sr.txt.gz')
                for sample in samples
            ],
            'PE_files': [
                str(sample.make_sv_evidence_path / f'{sample.id}.pe.txt.gz')
                for sample in samples
            ],
            'SD_files': [
                str(sample.make_sv_evidence_path / f'{sample.id}.sd.txt.gz')
                for sample in samples
            ],
            'ref_copy_number_autosomal_contigs': 2,
            'allosomal_contigs': ['chrX', 'chrY'],
            'gcnv_qs_cutoff': 30,
            'min_svsize': 50,
            'run_matrix_qc': True,
            'matrix_qc_distance': 1000000,
            'ref_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(sample.make_sv_evidence_path / f'{sample.id}.{caller}.vcf.gz')
                for sample in samples
            ]

        input_dict |= get_references(
            [
                'genome_file',
                'primary_contigs_fai',
                {'sd_locs_vcf': 'dbsnp_vcf'},
                {'cnmops_chrom_file': 'autosome_file'},
                'cnmops_exclude_list',
                {'cnmops_allo_file': 'allosome_file'},
                'cytoband',
                'mei_bed',
            ]
        )

        # reference panel gCNV models
        input_dict |= get_ref_panel()

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
                'cnmops_docker',
            ]
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=GatherBatchEvidence)
class ClusterBatch(CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#clusterbatch
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        * Clustered SV VCFs
        * Clustered depth-only call VCF
        * Metrics
        """

        ending_by_key = {
            'metrics_file_clusterbatch': 'metrics.tsv',
        }
        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'clustered_{caller}_vcf'] = f'clustered-{caller}.vcf.gz'
            ending_by_key[
                f'clustered_{caller}_vcf_index'
            ] = f'clustered-{caller}.vcf.gz.tbi'
        return {key: self.prefix / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """
        batch_evidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'del_bed': str(batch_evidence_d['merged_dels']),
            'dup_bed': str(batch_evidence_d['merged_dups']),
            'ped_file': str(make_combined_ped(cohort, self.prefix)),
            'depth_exclude_overlap_fraction': 0.5,
            'depth_interval_overlap': 0.8,
            'depth_clustering_algorithm': 'SINGLE_LINKAGE',
            'pesr_interval_overlap': 0.1,
            'pesr_breakend_window': 300,
            'pesr_clustering_algorithm': 'SINGLE_LINKAGE',
            'reference_fasta': str(get_fasta()),
            'reference_fasta_fai': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcf_tar'] = str(
                batch_evidence_d[f'std_{caller}_vcf_tar']
            )

        input_dict |= get_images(
            ['sv_base_mini_docker', 'sv_pipeline_docker', 'gatk_docker', 'linux_docker']
        )

        input_dict |= get_references(
            [
                {'contig_list': 'primary_contigs_list'},
                {'depth_exclude_intervals': 'depth_exclude_list'},
                {'pesr_exclude_intervals': 'pesr_exclude_list'},
            ]
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(CohortStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Metrics files
        """

        return {
            'metrics': self.prefix / 'metrics.tsv',
            'metrics_common': self.prefix / 'metrics_common.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)
        gatherbatchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'baf_metrics': gatherbatchevidence_d['merged_BAF'],
            'discfile': gatherbatchevidence_d['merged_PE'],
            'coveragefile': gatherbatchevidence_d['merged_bincov'],
            'splitfile': gatherbatchevidence_d['merged_SR'],
            'medianfile': gatherbatchevidence_d['median_cov'],
            'BAF_split_size': 10000,
            'RD_split_size': 10000,
            'PE_split_size': 10000,
            'SR_split_size': 1000,
            'common_cnv_size_cutoff': 5000,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'ref_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'sv_base_docker',
                'linux_docker',
            ]
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
                'rmsk',
                'segdups',
                {'autosome_contigs': 'autosome_file'},
                {'allosome_contigs': 'allosome_file'},
            ]
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[GenerateBatchMetrics, ClusterBatch])
class FilterBatch(CohortStage):
    """
    Filters poor quality variants and filters outlier samples. This workflow can
    be run all at once with the WDL at wdl/FilterBatch.wdl, or it can be run in three
    steps to enable tuning of outlier filtration cutoffs. The three sub-workflows are:

    * FilterBatchSites: Per-batch variant filtration
    * PlotSVCountsPerSample: Visualize SV counts per sample per type to help choose an
      IQR cutoff for outlier filtering, and preview outlier samples for a given cutoff
    * FilterBatchSamples: Per-batch outlier sample filtration; provide an appropriate
      outlier_cutoff_nIQR based on the SV count plots and outlier previews from step 2.
      Note that not removing high outliers can result in increased compute cost and
      a higher false positive rate in later steps.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        * Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        * Filtered depth-only call VCF with outlier samples excluded
        * Random forest cutoffs file
        * PED file with outlier samples excluded
        """

        ending_by_key: dict = {
            'metrics_file_filterbatch': 'metrics.tsv',
            'filtered_pesr_vcf': 'filtered_pesr_merged.vcf.gz',
            'cutoffs': 'cutoffs',
            'scores': 'updated_scores',
            'RF_intermediate_files': 'RF_intermediate_files.tar.gz',
            'outlier_samples_excluded_file': 'outliers.samples.list',
            'batch_samples_postOutlierExclusion_file': 'outliers_excluded.samples.list',
        }
        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'filtered_{caller}_vcf'] = f'filtered-{caller}.vcf.gz'

            # unsure why, scramble doesn't export this file
            if caller != 'scramble':
                ending_by_key[
                    f'sites_filtered_{caller}_vcf'
                ] = f'sites-filtered-{caller}.vcf.gz'

        ending_by_key['sv_counts'] = [
            f'{caller}.with_evidence.svcounts.txt' for caller in SV_CALLERS + ['depth']
        ]
        ending_by_key['sv_count_plots'] = [
            f'{caller}.with_evidence.all_SVTYPEs.counts_per_sample.png'
            for caller in SV_CALLERS + ['depth']
        ]
        d: dict[str, Path | list[Path]] = {}
        for key, ending in ending_by_key.items():
            if isinstance(ending, str):
                d[key] = self.prefix / ending
            elif isinstance(ending, list):
                d[key] = [self.prefix / e for e in ending]
        return d

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        metrics_d = inputs.as_dict(cohort, GenerateBatchMetrics)
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'evidence_metrics': metrics_d['metrics'],
            'evidence_metrics_common': metrics_d['metrics_common'],
            'outlier_cutoff_nIQR': '6',
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'linux_docker',
            ]
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
            ]
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage
class MergeBatchSites(CohortStage):
    """
    I'm adding these here as it logically sits between FilterBatch and GenotypeBatch
    Though it runs independently, with input and output for this stage communicated
    via config (at least until the more complex batching/cohort logic comes along)
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        generate them there outputs
        """
        return {
            'cohort_pesr_vcf': self.prefix / 'cohort_pesr.vcf.gz',
            'cohort_depth_vcf': self.prefix / 'cohort_depth.vcf.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        generate a MergeBatchSites job
        """

        batch_names = get_config()['workflow']['batch_names']
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        pesr_vcfs = [
            batch_prefix / batch_name / 'FilterBatch' / 'filtered_pesr_merged.vcf.gz'
            for batch_name in batch_names
        ]
        depth_vcfs = [
            batch_prefix / batch_name / 'FilterBatch' / 'filtered-depth.vcf.gz'
            for batch_name in batch_names
        ]

        input_dict: dict[str, Any] = {
            'cohort': cohort.name,
            'depth_vcfs': depth_vcfs,
            'pesr_vcfs': pesr_vcfs,
        }
        input_dict |= get_images(['sv_pipeline_docker'])
        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[FilterBatch, GatherBatchEvidence])
class GenotypeBatch(CohortStage):
    """
    Genotypes a batch of samples across filtered variants combined across all batches.
    This requires a separate stage shoved in between FilterBatch and GenotypeBatch
    to aggregate all the filterBatch outputs across all sub-batches
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """
        ending_by_key = {
            'sr_bothside_pass': 'genotype_SR_part2_bothside_pass.txt',
            'sr_background_fail': 'genotype_SR_part2_background_fail.txt',
            'trained_PE_metrics': 'pe_metric_file.txt',
            'trained_SR_metrics': 'sr_metric_file.txt',
            'regeno_coverage_medians': 'regeno.coverage_medians_merged.bed',
            'metrics_file_genotypebatch': 'metrics.tsv',
        }

        for mode in ['pesr', 'depth']:
            ending_by_key |= {
                f'trained_genotype_{mode}_pesr_sepcutoff': f'{mode}.pesr_sepcutoff.txt',
                f'trained_genotype_{mode}_depth_sepcutoff': f'{mode}.depth_sepcutoff.txt',
                f'genotyped_{mode}_vcf': f'{mode}.vcf.gz',
                f'genotyped_{mode}_vcf_index': f'{mode}.vcf.gz.tbi',
            }

        return {key: self.prefix / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        filterbatch_d = inputs.as_dict(cohort, FilterBatch)
        batchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'n_per_split': 5000,
            'n_RD_genotype_bins': 100000,
            'coveragefile': batchevidence_d['merged_bincov'],  # unsure
            'coveragefile_index': batchevidence_d['merged_bincov_index'],  # unsure
            'discfile': batchevidence_d['merged_PE'],
            'discfile_index': batchevidence_d['merged_PE_index'],
            'splitfile': batchevidence_d['merged_SR'],
            'splitfile_index': batchevidence_d['merged_SR_index'],
            'medianfile': batchevidence_d['median_cov'],
            'rf_cutoffs': filterbatch_d['cutoffs'],
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'reference_build': 'hg38',
        }

        # pull out the merged VCF from MergeBatchSites
        for mode in ['pesr', 'depth']:
            input_dict[f'batch_{mode}_vcf'] = filterbatch_d[f'filtered_{mode}_vcf']
            input_dict[f'cohort_{mode}_vcf'] = get_config()['workflow'][
                f'cohort_{mode}_vcf'
            ]

        input_dict |= get_images(
            ['sv_pipeline_docker', 'sv_base_mini_docker', 'linux_docker']
        )
        input_dict |= get_references(
            ['primary_contigs_list', 'bin_exclude', 'seed_cutoffs', 'pesr_exclude_list']
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)
