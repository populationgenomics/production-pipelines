"""
The first multi-sample workflow, containing stages from
GenerateBatchMetrics to Genotype Batch

Though, there is a separate process required to merge files
across batches between FilterBatch and GenotypeBatch
"""

from typing import Any

from cpg_utils import Path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, try_get_ar_guid
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    SV_CALLERS,
    CromwellJobSizes,
    _sv_batch_meta,
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_ref_panel,
    get_references,
    make_combined_ped,
)
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)


@stage
class GatherBatchEvidence(CohortStage):
    """
    This is the first Stage in the multisample GATK-SV workflow, running on a
    controlled subset of SGs as determined by the output of CreateSampleBatches.

    Using Analysis-Runner, include three additional config files:

    - configs/gatk_sv/stop_at_filter_batch.toml
        - this contains the instruction to stop prior to FilterBatch
    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - A custom config with the specific cohort in workflow.only_sgs

    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl

    it's critical to separate the ending with a dot, e.g.: `*.sr.txt.gz`,
    These files are passed to `gatk PrintSVEvidence`, that determines file
    format based on the file name.
    It would strongly expect the files to end exactly with either
    `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
    "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
    no suitable codecs found".
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
            'merged_SR': f'{self.name}.sr.txt.gz',
            'merged_SR_index': f'{self.name}.sr.txt.gz.tbi',
            'merged_PE': f'{self.name}.pe.txt.gz',
            'merged_PE_index': f'{self.name}.pe.txt.gz.tbi',
            'merged_BAF': f'{self.name}.baf.txt.gz',
            'merged_BAF_index': f'{self.name}.baf.txt.gz.tbi',
            'merged_bincov': f'{self.name}.RD.txt.gz',
            'merged_bincov_index': f'{self.name}.RD.txt.gz.tbi',
            'median_cov': 'medianCov.transposed.bed',
            'merged_dels': 'DEL.bed.gz',
            'merged_dups': 'DUP.bed.gz',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_matrix_qc'], False):
            ending_by_key.update(
                {
                    'Matrix_QC_plot': '00_matrix_FC_QC.png',
                    'SR_stats': 'SR.QC_matrix.txt',
                    'PE_stats': 'PE.QC_matrix.txt',
                    'BAF_stats': 'BAF.QC_matrix.txt',
                    'RD_stats': 'RD.QC_matrix.txt',
                },
            )

        for caller in SV_CALLERS:
            ending_by_key[f'std_{caller}_vcf_tar'] = f'{caller}.tar.gz'

        return {key: self.prefix / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """Add jobs to Batch"""
        sequencing_groups = cohort.get_sequencing_groups(only_active=True)

        input_dict: dict[str, Any] = {
            'batch': get_workflow().output_version,
            'samples': [sg.id for sg in sequencing_groups],
            'ped_file': str(make_combined_ped(cohort, self.prefix)),
            'counts': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.coverage_counts.tsv.gz')
                for sequencing_group in sequencing_groups
            ],
            'SR_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sr.txt.gz')
                for sequencing_group in sequencing_groups
            ],
            'PE_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.pe.txt.gz')
                for sequencing_group in sequencing_groups
            ],
            'SD_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sd.txt.gz')
                for sequencing_group in sequencing_groups
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
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.{caller}.vcf.gz')
                for sequencing_group in sequencing_groups
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
            ],
        )

        # reference panel gCNV models
        input_dict |= get_ref_panel()

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
                'cnmops_docker',
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # this step runs for approximately 15 hours
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.LARGE,
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

        ending_by_key = {}

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_module_metrics'], False):
            ending_by_key['metrics_file_clusterbatch'] = 'metrics.tsv'

        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'clustered_{caller}_vcf'] = f'clustered-{caller}.vcf.gz'
            ending_by_key[f'clustered_{caller}_vcf_index'] = f'clustered-{caller}.vcf.gz.tbi'
        return {key: self.prefix / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """
        batch_evidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': get_workflow().output_version,
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
            input_dict[f'{caller}_vcf_tar'] = str(batch_evidence_d[f'std_{caller}_vcf_tar'])

        input_dict |= get_images(
            [
                'linux_docker',
                'sv_pipeline_base_docker',
                'gatk_docker',
                'sv_base_mini_docker',
                'sv_pipeline_docker',
            ],
        )

        input_dict |= get_references(
            [
                {'contig_list': 'primary_contigs_list'},
                {'depth_exclude_intervals': 'depth_exclude_list'},
                {'pesr_exclude_intervals': 'pesr_exclude_list'},
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for approx 1 hour
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(CohortStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, str | Path]:
        """
        Metrics files
        """

        outputs = {
            'metrics': self.prefix / 'metrics.tsv',
            'metrics_common': self.prefix / 'metrics_common.tsv',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_module_metrics'], False):
            outputs['metrics_file_batchmetrics'] = f'GenerateBatchMetrics.{get_workflow().output_version}'

        return outputs

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)
        gatherbatchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': get_workflow().output_version,
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
                'sv_pipeline_rdtest_docker',
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_base_docker',
                'linux_docker',
            ],
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
                'rmsk',
                'segdups',
                {'autosome_contigs': 'autosome_file'},
                {'allosome_contigs': 'allosome_file'},
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for approx 4-5 hours
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[GenerateBatchMetrics, ClusterBatch])
class FilterBatch(CohortStage):
    """
    Filters poor quality variants and filters outlier samples.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        * Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        * Filtered depth-only call VCF with outlier samples excluded
        * Random forest cutoffs file
        * PED file with outlier samples excluded
        """

        ending_by_key: dict = {
            'filtered_pesr_vcf': 'filtered_pesr_merged.vcf.gz',
            'cutoffs': 'cutoffs',
            'scores': 'updated_scores',
            'RF_intermediate_files': 'RF_intermediate_files.tar.gz',
            'outlier_samples_excluded_file': 'outliers.samples.list',
            'batch_samples_postOutlierExclusion_file': 'outliers_excluded.samples.list',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_module_metrics'], False):
            ending_by_key['metrics_file_filterbatch'] = 'metrics.tsv'

        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'filtered_{caller}_vcf'] = f'filtered-{caller}.vcf.gz'

            # unsure why, scramble doesn't export this file
            if caller != 'scramble':
                ending_by_key[f'sites_filtered_{caller}_vcf'] = f'sites-filtered-{caller}.vcf.gz'

        ending_by_key['sv_counts'] = [f'{caller}.with_evidence.svcounts.txt' for caller in SV_CALLERS + ['depth']]
        ending_by_key['sv_count_plots'] = [
            f'{caller}.with_evidence.all_SVTYPEs.counts_per_sample.png' for caller in SV_CALLERS + ['depth']
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
            'batch': get_workflow().output_version,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'evidence_metrics': metrics_d['metrics'],
            'evidence_metrics_common': metrics_d['metrics_common'],
            'outlier_cutoff_nIQR': '6',
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            ['sv_pipeline_base_docker', 'sv_pipeline_docker', 'sv_base_mini_docker', 'linux_docker'],
        )

        input_dict |= get_references(['primary_contigs_list'])

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for over an hour
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage
class MergeBatchSites(CohortStage):
    """
    This Stage runs separately between GenerateBatchMetrics and FilterBatch
    This takes the component VCFs from individual batches and merges them into
    a single VCF (one for PESR and one for depth-only calls).

    The single-stage workflow used here is gatk_sv_sandwich.

    This only has to be run once for all sub-cohorts, and technically doesn't
    need to be associated with a group of SGs. However, including one (via the
    config containing `only_sgs`) allows us to keep the output VCFs within the
    output structure generated by previous stages in this run.

    Using Analysis-Runner, include three additional config files:
    - configs/gatk_sv/all_batch_names.toml
        - add all sub-cohort hashes to the batch_names list
        - for ongoing future runs, we may want to include ALL prior hashes
    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - A custom config with a specific cohort in workflow.only_sgs
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

        batch_names = config_retrieve(['workflow', 'batch_names'])
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        pesr_vcfs = [
            batch_prefix / batch_name / 'FilterBatch' / 'filtered_pesr_merged.vcf.gz' for batch_name in batch_names
        ]
        depth_vcfs = [batch_prefix / batch_name / 'FilterBatch' / 'filtered-depth.vcf.gz' for batch_name in batch_names]

        input_dict: dict[str, Any] = {'cohort': cohort.name, 'depth_vcfs': depth_vcfs, 'pesr_vcfs': pesr_vcfs}
        input_dict |= get_images(['sv_pipeline_docker'])
        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(
    required_stages=[FilterBatch, GatherBatchEvidence],
    analysis_type='sv',
    analysis_keys=[f'genotyped_{mode}_vcf' for mode in ['pesr', 'depth']],
    update_analysis_meta=_sv_batch_meta,
)
class GenotypeBatch(CohortStage):
    """
    Genotypes a batch of samples across filtered variants combined across all batches.
    In the current WF setup, this Stage is run only once MergeBatchSites (via the WF
    `gatk_sv_sandwich`) has completed, and a pair of VCFs containing all merged
    FilterBatch outputs individual batches has been generated.

    Using Analysis-Runner, include three additional config files:

    - configs/gatk_sv/genotypebatch.toml
        - this contains the instruction to only run GenotypeBatch
        - update the cohort_depth_vcf and cohort_pesr_vcf values
    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - A custom config with the specific cohort in workflow.only_sgs
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """

        # this workflow requires re-genotyping previous batches with each new
        # expansion of the cohort.
        # we want to keep all versions of these outputs generated from separate
        # inputs, so we incorporate the hash from the MergeBatchSites file into
        # the files generated by this step to differentiate the outputs
        if mbs_file := config_retrieve(['workflow', 'cohort_depth_vcf'], False):
            cohort_hash = mbs_file.split('/')[4][-10:]
        else:
            # this method is still run, even if the Stage is skipped
            # adding a failure here prevents the need for a dummy config variable
            print('cohort_depth_vcf not found in config, wont run GenotypeBatch')
            return {}

        ending_by_key = {
            'sr_bothside_pass': 'genotype_SR_part2_bothside_pass.txt',
            'sr_background_fail': 'genotype_SR_part2_background_fail.txt',
            'trained_PE_metrics': 'pe_metric_file.txt',
            'trained_SR_metrics': 'sr_metric_file.txt',
            'regeno_coverage_medians': 'regeno.coverage_medians_merged.bed',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_module_metrics'], False):
            ending_by_key['metrics_file_genotypebatch'] = 'metrics.tsv'

        for mode in ['pesr', 'depth']:
            ending_by_key |= {
                f'trained_genotype_{mode}_pesr_sepcutoff': f'{mode}.pesr_sepcutoff.txt',
                f'trained_genotype_{mode}_depth_sepcutoff': f'{mode}.depth_sepcutoff.txt',
                f'genotyped_{mode}_vcf': f'{mode}.vcf.gz',
                f'genotyped_{mode}_vcf_index': f'{mode}.vcf.gz.tbi',
            }

        return {key: self.prefix / f'{cohort_hash}_{fname}' for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        # if this is run before the Sandwich Stage, we can't run GenotypeBatch
        if not config_retrieve(['workflow', 'cohort_depth_vcf'], False):
            print('cohort_depth_vcf not found in config, wont run GenotypeBatch')
            return self.make_outputs(cohort, skipped=True)

        filterbatch_d = inputs.as_dict(cohort, FilterBatch)
        batchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': get_workflow().output_version,
            'n_per_split': 5000,
            'n_RD_genotype_bins': 100000,
            'coveragefile': batchevidence_d['merged_bincov'],
            'coveragefile_index': batchevidence_d['merged_bincov_index'],
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
            input_dict[f'cohort_{mode}_vcf'] = config_retrieve(['workflow', f'cohort_{mode}_vcf'])

        input_dict |= get_images(
            [
                'sv_pipeline_base_docker',
                'sv_base_mini_docker',
                'sv_pipeline_docker',
                'sv_pipeline_rdtest_docker',
                'linux_docker',
            ],
        )
        input_dict |= get_references(['primary_contigs_list', 'bin_exclude', 'seed_cutoffs', 'pesr_exclude_list'])
        # 2 additional (optional) references CPG doesn't have:
        # - sr_hom_cutoff_multiplier
        # - sr_median_hom_ins

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)
