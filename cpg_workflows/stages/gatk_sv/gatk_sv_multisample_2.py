"""
The second multi-sample workflow, containing the following stages:
MakeCohortVCF and AnnotateVCF
"""


from typing import Any

from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import stage, StageOutput, StageInput, Cohort, CohortStage

from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
    make_combined_ped,
    _sv_batch_meta,
    _sv_filtered_meta,
    SV_CALLERS,
)
from cpg_workflows.jobs import ploidy_table_from_ped


@stage
class MakeCohortVcf(CohortStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    This Stage is the first in the gatk_sv_multisample_2 workflow. Using Analysis-Runner,
    include two additional config files:

    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - configs/gatk_sv/all_batch_names.toml
        - add all sub-cohort hashes to the batch_names list
    - A custom config with all combined SGs in workflow.only_sgs

    Once this becomes a more frequently run process, we'll have to investigate whether
    we want to include ALL prior batches, or just the batches which were run in the
    earlier pipeline steps this time around. That may depend on how much these stages
    cost, and whether generating the final cross-batch VCF/ES-index is possible from
    several component joint calls, or we require one single joint call across all SGs
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """create output paths"""
        return {
            'vcf': self.prefix / 'cleaned.vcf.gz',
            'vcf_index': self.prefix / 'cleaned.vcf.gz.tbi',
            'vcf_qc': self.prefix / 'cleaned_SV_VCF_QC_output.tar.gz',
            'metrics_file_makecohortvcf': self.prefix / 'metrics.tsv',
            # if merge_intermediate_vcfs is enabled
            # 'cluster_vcf': '.combine_batches.vcf.gz',
            # 'cluster_vcf_index': '.combine_batches.vcf.gz.tbi',
            # 'complex_resolve_vcf': '.complex_resolve.vcf.gz',
            # 'complex_resolve_vcf_index': '.complex_resolve.vcf.gz.tbi',
            # 'complex_genotype_vcf': '.complex_genotype.vcf.gz',
            # 'complex_genotype_vcf_index': '.complex_genotype.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        This is a little bit spicy. Instead of taking a direct dependency on the
        previous stage, we use the output hash to find all the previous batches

        Replacing this nasty mess with nested/fancy cohorts would be ace
        """

        batch_names = get_config()['workflow']['batch_names']
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        pesr_vcfs = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'pesr.vcf.gz'
            for batch_name in batch_names
        ]
        depth_vcfs = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'depth.vcf.gz'
            for batch_name in batch_names
        ]
        sr_pass = [
            batch_prefix
            / batch_name
            / 'GenotypeBatch'
            / 'genotype_SR_part2_bothside_pass.txt'
            for batch_name in batch_names
        ]
        sr_fail = [
            batch_prefix
            / batch_name
            / 'GenotypeBatch'
            / 'genotype_SR_part2_background_fail.txt'
            for batch_name in batch_names
        ]
        depth_depth_cutoff = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'depth.depth_sepcutoff.txt'
            for batch_name in batch_names
        ]
        filter_batch_cutoffs = [
            batch_prefix / batch_name / 'FilterBatch' / 'cutoffs'
            for batch_name in batch_names
        ]
        bincov_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'GatherBatchEvidence.RD.txt.gz'
            for batch_name in batch_names
        ]
        disc_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'GatherBatchEvidence.pe.txt.gz'
            for batch_name in batch_names
        ]
        median_cov_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'medianCov.transposed.bed'
            for batch_name in batch_names
        ]

        input_dict: dict[str, Any] = {
            'cohort_name': cohort.name,
            'batches': batch_names,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'chr_x': 'chrX',
            'chr_y': 'chrY',
            'min_sr_background_fail_batches': 0.5,
            'max_shard_size_resolve': 500,
            'max_shards_per_chrom_clean_vcf_step1': 200,
            'min_records_per_shard_clean_vcf_step1': 5000,
            'clean_vcf1b_records_per_shard': 10000,
            'samples_per_clean_vcf_step2_shard': 100,
            'clean_vcf5_records_per_shard': 5000,
            'random_seed': 0,
            # not explicit, but these VCFs require indices
            'pesr_vcfs': pesr_vcfs,
            'depth_vcfs': depth_vcfs,
            'disc_files': disc_files,
            'bincov_files': bincov_files,
            'raw_sr_bothside_pass_files': sr_pass,
            'raw_sr_background_fail_files': sr_fail,
            'depth_gt_rd_sep_files': depth_depth_cutoff,
            'median_coverage_files': median_cov_files,
            'rf_cutoff_files': filter_batch_cutoffs,
        }

        input_dict |= get_references(
            [
                'bin_exclude',
                'mei_bed',
                'depth_exclude_list',
                'empty_file',
                # same attr, two names
                'primary_contigs_list',
                {'contig_list': 'primary_contigs_list'},
                {'allosome_fai': 'allosome_file'},
                {'cytobands': 'cytoband'},
                {'pe_exclude_list': 'pesr_exclude_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_pipeline_base_docker',
                'sv_pipeline_hail_docker',
                'sv_pipeline_updates_docker',
                'sv_pipeline_rdtest_docker',
                'sv_pipeline_qc_docker',
                'sv_base_mini_docker',
                'linux_docker',
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


@stage(required_stages=MakeCohortVcf)
class FormatVcfForGatk(CohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """


        Args:
            cohort (Cohort): cohort of all samples (across several sub-cohort batches)
            inputs (StageInput): access to prior inputs
        """

        make_vcf_d = inputs.as_dict(cohort, MakeCohortVcf)
        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'vcf': make_vcf_d['vcf'],
            'ped_file': make_combined_ped(cohort, self.prefix),
        }
        input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=MakeCohortVcf)
class JoinRawCalls(CohortStage):
    """
    Joins all individually clustered caller results
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'joined_raw_calls_vcf': self.prefix / 'raw_clustered_calls.vcf.gz',
            'joined_raw_calls_vcf_index': self.prefix
            / 'raw_clustered_calls.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """ """

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'reference_fasta': get_fasta(),
            'reference_fasta_fai': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
        }
        input_dict |= get_images(
            ['gatk_docker', 'sv_pipeline_docker', 'sv_base_mini_docker']
        )
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        # add all clustered _caller_ files, plus indices
        batch_names = get_config()['workflow']['batch_names']
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        for caller in SV_CALLERS + ['depth']:
            input_dict[f'clustered_{caller}_vcfs'] = [
                batch_prefix
                / batch_name
                / 'ClusterBatch'
                / f'clustered-{caller}.vcf.gz'
                for batch_name in batch_names
            ]
            input_dict[f'clustered_{caller}_vcf_indexes'] = [
                batch_prefix
                / batch_name
                / 'ClusterBatch'
                / f'clustered-{caller}.vcf.gz.tbi'
                for batch_name in batch_names
            ]
        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[JoinRawCalls, FormatVcfForGatk])
class SVConcordance(CohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'concordance_vcf': self.prefix / 'sv_concordance.vcf.gz',
            'concordance_vcf_index': self.prefix / 'sv_concordance.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV concordance
        """
        raw_calls = inputs.as_dict(cohort, JoinRawCalls)
        format_vcf = inputs.as_dict(cohort, FormatVcfForGatk)

        input_dict: dict[str, Any] = {
            'output_prefix': cohort.name,
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'eval_vcf': format_vcf['gatk_formatted_vcf'],
            'truth_vcf': raw_calls['joined_raw_calls_vcf'],
        }
        input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=SVConcordance)
class GeneratePloidyTable(CohortStage):
    """
    Quick PythonJob to generate a ploidy table
    Calls a homebrewed version of this table generator:
    github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/ploidy_table_from_ped.py
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        only one output, the ploidy table
        """

        return {'ploidy_table': self.prefix / 'ploidy_table.txt'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        py_job = get_batch().new_python_job('create_ploidy_table')
        py_job.image(get_config()['workflow']['driver_image'])

        ped_path = make_combined_ped(cohort, self.prefix)
        contig_path = get_references(['primary_contigs_list'])['primary_contigs_list']

        expected_d = self.expected_outputs(cohort)
        py_job.call(
            ploidy_table_from_ped.generate_ploidy_table,
            ped_path,
            contig_path,
            expected_d['ploidy_table'],
        )

        return self.make_outputs(cohort, data=expected_d, jobs=py_job)


@stage(
    required_stages=[GeneratePloidyTable, SVConcordance],
    analysis_type='sv',
    analysis_keys=['filtered_vcf'],
    update_analysis_meta=_sv_filtered_meta,
)
class FilterGenotypes(CohortStage):
    """
    Steps required to post-filter called genotypes
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'filtered_vcf': self.prefix / 'filtered.vcf.gz',
            'filtered_vcf_index': self.prefix / 'filtered.vcf.gz.tbi',
            'main_vcf_qc_tarball': self.prefix / 'filtered_SV_VCF_QC_output.tar.gz',
            'unfiltered_recalibrated_vcf': self.prefix
            / 'unfiltered_recalibrated.vcf.gz',
            'unfiltered_recalibrated_vcf_index': self.prefix
            / 'unfiltered_recalibrated.vcf.gz.tbi',
            'vcf_optimization_table': self.prefix / 'vcf_optimization_table.tsv.gz',
            'sl_cutoff_qc_tarball': self.prefix / 'sl_cutoff_SV_VCF_QC_output.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        input_dict = {
            'output_prefix': cohort.name,
            'vcf': inputs.as_dict(cohort, SVConcordance)['concordance_vcf'],
            'ploidy_table': inputs.as_dict(cohort, GeneratePloidyTable)['ploidy_table'],
        }

        input_dict |= get_images(
            ['gatk_docker', 'linux_docker', 'sv_base_mini_docker', 'sv_pipeline_docker']
        )

        input_dict |= get_references(
            [{'gq_recalibrator_model_file': 'aou_filtering_model'}]
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


@stage(
    required_stages=FilterGenotypes,
    analysis_type='sv',
    analysis_keys=['output_vcf'],
    update_analysis_meta=_sv_batch_meta,
)
class AnnotateVcf(CohortStage):
    """
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        return {
            'output_vcf': self.prefix / 'filtered_annotated.vcf.gz',
            'output_vcf_idx': self.prefix / 'filtered_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        make_vcf_d = inputs.as_dict(cohort, FilterGenotypes)

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'vcf': make_vcf_d['filtered_vcf'],
            'vcf_idx': make_vcf_d['filtered_vcf_index'],
            'ped_file': make_combined_ped(cohort, self.prefix),
            'sv_per_shard': 5000,
            'max_shards_per_chrom_step1': 200,
            'min_records_per_shard_step1': 5000,
        }
        input_dict |= get_references(
            [
                'protein_coding_gtf',
                {'contig_list': 'primary_contigs_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'gatk_docker',
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
