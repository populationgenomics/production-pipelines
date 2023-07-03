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
)


@stage
class MakeCohortVcf(CohortStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    If RegenotypeCNVs is run, this stage will use the output of that stage as input.
    Otherwise, it will use the output of GenotypeBatch.
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
            batch_prefix / batch_name / 'GatherBatchEvidence' / 'GatherBatchEvidence.RD.txt.gz'
            for batch_name in batch_names
        ]
        disc_files = [
            batch_prefix / batch_name / 'GatherBatchEvidence' / 'GatherBatchEvidence.pe.txt.gz'
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
            'batches': [batch_names],
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


@stage(
    required_stages=MakeCohortVcf,
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
        # vcf_path = str(self.tmp_prefix / 'sv' / f'{dataset.name}.vcf.gz')
        return {
            'output_vcf': self.prefix / 'annotated.vcf.gz',
            'output_vcf_idx': self.prefix / 'annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        make_vcf_d = inputs.as_dict(cohort, MakeCohortVcf)

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'vcf': make_vcf_d['vcf'],
            'vcf_idx': make_vcf_d['vcf_index'],
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
