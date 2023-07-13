"""
Stage to call SNVs in the mitochondrial genome of a single sequencing_group.

Reimplemented version of;
https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

"""
from typing import Any

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import mito, picard, mito_cohort, joint_genotyping, vep, seqr_loader
from cpg_workflows.stages.align import Align
from cpg_workflows.stages.cram_qc import CramQC
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import exists
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroupStage,
    SequencingGroup,
    CohortStage,
)


MITO_REF = {
    'dict': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict',
    'base': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta',
    'amb': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb',
    'ann': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann',
    'bwt': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt',
    'fai': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai',
    'pac': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac',
    'sa': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa',
}

SHIFTED_MITO_REF = {
    'dict': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict',
    'base': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta',
    'amb': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb',
    'ann': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann',
    'bwt': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt',
    'fai': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai',
    'pac': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac',
    'sa': 'gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa',
    'shift_back_chain': 'gs://cpg-common-main/references/hg38/v0/chrM/ShiftBack.chain',
}

CONTROL_REGION_INTERVALS = {
    'control_region_shifted': 'gs://cpg-common-main/references/hg38/v0/chrM/control_region_shifted.chrM.interval_list',
    'non_control_region': 'gs://cpg-common-main/references/hg38/v0/chrM/non_control_region.chrM.interval_list',
}


@stage(
    required_stages=Align,
    # TODO: add suitable analysis types
)
class RealignMito(SequencingGroupStage):
    """
    Re-align mitochondrial genome of a single sequencing_group.

    This is a re-implementation of the broad pipeline (as of 03/22) used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

    Mitochondrial variant calling is a subtle art with potential for artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to avoid these
    issues, this stage has blind faith in the Broad pipeline and faithfully re-implements
    as much of the logic, tools and configuration as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sequencing_group cram.
        - Realignment of chrM reads to chrM reference using bwa.
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted"
            reference is the same sequence but with a different linearisation point. This
            is used to overcome mapping artifacts at the boundary (aka the "control
            region" of the mt genome).
        - Generation of overall and  base level alignment statistics.

    Mitochondrial reference indexes and lit over chains were sourced from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).

    Requires:
        sequencing_group cram from Align stage.

    Outputs:
        non_shifted_cram: Sorted cram with duplicates marked aligned to the standard
            hg38 chrM.
        shifted_cram: Sorted cram with duplicates marked aligned to a copy of chrM
            with the linearisation location "shifted"
        base_level_coverage_metrics: per base coverage needed to differentiate between
            0/0 and ./. genotypes until mutect can generate a gVCF.
        coverage_metrics: CollectWgsMetrics output.
        theoretical_sensitivity_metrics: CollectWgsMetrics output.

    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        analysis = sequencing_group.dataset.analysis_prefix()
        return {
            'non_shifted_cram': main / 'mito' / f'{sequencing_group.id}.mito.cram',
            'shifted_cram': main / 'mito' / f'{sequencing_group.id}.shifted_mito.cram',
            'base_level_coverage_metrics': main
            / 'mito'
            / f'{sequencing_group.id}.base_level_coverage.tsv',
            'coverage_metrics': analysis
            / 'mito'
            / f'{sequencing_group.id}.coverage_metrics.txt',
            'theoretical_sensitivity_metrics': analysis
            / 'mito'
            / f'{sequencing_group.id}.theoretical_sensitivity.txt',
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # Mitochondrial specific reference files.
        mito_ref = get_batch().read_input_group(**MITO_REF)
        shifted_mito_ref = get_batch().read_input_group(**SHIFTED_MITO_REF)
        intervals = get_batch().read_input_group(**CONTROL_REGION_INTERVALS)

        jobs = []

        # Extract reads mapped to chrM
        cram_path = inputs.as_path(sequencing_group, Align, 'cram')
        crai_path = inputs.as_path(sequencing_group, Align, 'crai')
        subset_bam_j = mito.subset_cram_to_chrM(
            b=get_batch(),
            cram_path=CramPath(cram_path, crai_path),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(subset_bam_j)
        assert isinstance(subset_bam_j.output_bam, hb.ResourceGroup)

        # Align extracted reads to chrM using bwa
        realign_j = mito.mito_realign(
            b=get_batch(),
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(realign_j)
        assert isinstance(realign_j.output_cram, hb.ResourceFile)

        # Mark duplicates (non-shifted)
        realign_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=realign_j.output_cram,
            fasta_reference=mito_ref,
            output_path=self.expected_outputs(sequencing_group)['non_shifted_cram'],
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=True,
        )
        jobs.append(realign_mkdup_j)
        assert isinstance(realign_mkdup_j, Job)
        assert isinstance(realign_mkdup_j.output_cram, hb.ResourceGroup)

        # Align extracted reads to "shifted" chrM using bwa
        shifted_realign_j = mito.mito_realign(
            b=get_batch(),
            sequencing_group_id=sequencing_group.id,
            input_bam=subset_bam_j.output_bam,
            mito_ref=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_realign_j)
        assert isinstance(shifted_realign_j.output_cram, hb.ResourceFile)

        # Mark duplicates (shifted)
        shifted_mkdup_j = picard.markdup(
            b=get_batch(),
            sorted_bam=shifted_realign_j.output_cram,
            fasta_reference=shifted_mito_ref,
            output_path=self.expected_outputs(sequencing_group)['shifted_cram'],
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=True,
        )
        jobs.append(shifted_mkdup_j)
        assert isinstance(shifted_mkdup_j, Job)
        assert isinstance(shifted_mkdup_j.output_cram, hb.ResourceGroup)

        # Collect coverage metrics (only on non-shifted)
        coverage_metrics_J = mito.collect_coverage_metrics(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            reference=mito_ref,
            metrics=self.expected_outputs(sequencing_group)['coverage_metrics'],
            theoretical_sensitivity=self.expected_outputs(sequencing_group)[
                'theoretical_sensitivity_metrics'
            ],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(coverage_metrics_J)

        # CoverageAtEveryBase
        # This is a temporary task to handle "joint calling" until Mutect2 can produce a
        # GVCF. This proivdes coverage at each base so low coverage sites can be
        # considered ./. rather than 0/0.
        non_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            intervals_list=intervals.non_control_region,
            reference=mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(non_control_region_coverage_j)
        assert isinstance(
            non_control_region_coverage_j.per_base_coverage, hb.ResourceFile
        )

        shifted_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=shifted_mkdup_j.output_cram,
            intervals_list=intervals.control_region_shifted,
            reference=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_control_region_coverage_j)
        assert isinstance(
            shifted_control_region_coverage_j.per_base_coverage, hb.ResourceFile
        )

        # Merge coverage stats
        merge_coverage_j = mito.merge_coverage(
            b=get_batch(),
            non_cr_coverage=non_control_region_coverage_j.per_base_coverage,
            shifted_cr_coverage=shifted_control_region_coverage_j.per_base_coverage,
            merged_coverage=self.expected_outputs(sequencing_group)[
                'base_level_coverage_metrics'
            ],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_coverage_j)

        return self.make_outputs(
            sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs
        )


@stage(
    required_stages=[RealignMito, CramQC]
    # TODO: add suitable analysis types
)
class GenotypeMito(SequencingGroupStage):
    """
    Call SNVs in the mitochondrial genome of a single sequencing_group.
    This is a re-implementation of the broad pipeline as of 03/22 that was used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json
    Mitochondrial variant calling is a subtle art with potential for artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to avoid these
    issues, this stage has blind faith in the gnomAD pipeline and faithfully re-implements
    as much of the logic, tools and configuration as possible.
    The main phases of analysis include:
        - Calling of variants from non-shifted cram using mutect2
        - Calling of variants in control region in shifted cram using mutect2
        - Merging of the two call sets into a single vcf in normal chrM coordinate space.
        - Variant filtering part 1: Exclude black list of known problem sites and high
            alt allele counts.
        - Estimate contamination with haplocheckCLI (using known mito haplotypes)
        - Variant filtering part 2: exclude variants with VAF below contamination
            estimate.
        - Export final vcf
    Mitochondrial reference indexes and filtering black lists were copied from Broad
    (gs://gcp-public-data--broad-references/hg38/v0/chrM/).
    Requires:
        sequencing_group cram from Align stage.
    Outputs:
        out_vcf: the final filtered vcf for downstream use
        haplocheck_metrics: Metrics generated from the haplocheckCLI tool including an
            estimate of contamination and the predicted mitochondrial haplotype found.
    Configuration options:
    The following are surfaced as configurable parameters in the Broad WDL. Other
    parameters hardcoded in the WDL are also hardcoded in this pipeline.
        mito_snv.vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
        mito_snv.f_score_beta: "F-Score beta balances the filtering strategy between
            recall and precision. The relative weight of recall to precision."
    Not Implemented:
        - The Broad wdl allows for use of verifyBamID as a second input for contamination
            estimation. This has not been implemented yet but is probably a good idea.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        analysis = sequencing_group.dataset.analysis_prefix()
        return {
            'out_vcf': main / 'mito' / f'{sequencing_group.id}.mito.vcf.bgz',
            'haplocheck_metrics': analysis
            / 'mito'
            / f'{sequencing_group.id}.haplocheck.txt',
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # Mitochondrial specific reference files.
        mito_ref = get_batch().read_input_group(**MITO_REF)
        shifted_mito_ref = get_batch().read_input_group(**SHIFTED_MITO_REF)

        jobs = []

        # Get input resources
        non_shifted_cram = get_batch().read_input_group(
            cram=str(inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram')),
            crai=str(inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram'))
            + '.crai',
        )
        shifted_cram = get_batch().read_input_group(
            cram=str(inputs.as_path(sequencing_group, RealignMito, 'shifted_cram')),
            crai=str(inputs.as_path(sequencing_group, RealignMito, 'shifted_cram'))
            + '.crai',
        )
        if get_config()['mito_snv']['use_verifybamid']:
            verifybamid_output = get_batch().read_input(
                str(inputs.as_path(sequencing_group, CramQC, 'verify_bamid')),
            )
        else:
            verifybamid_output = None

        # Call variants on WT genome
        call_j = mito.mito_mutect2(
            b=get_batch(),
            cram=non_shifted_cram,
            reference=mito_ref,
            region='chrM:576-16024',  # Exclude the control region.
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(call_j)
        assert isinstance(call_j.output_vcf, hb.ResourceGroup)

        # Call variants in ONLY the control region using the shifted reference
        shifted_call_j = mito.mito_mutect2(
            b=get_batch(),
            cram=shifted_cram,
            reference=shifted_mito_ref,
            region='chrM:8025-9144',  # Only call inside the control region.
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_call_j)
        assert isinstance(shifted_call_j.output_vcf, hb.ResourceGroup)

        # Merge the wt and shifted VCFs
        merge_j = mito.liftover_and_combine_vcfs(
            b=get_batch(),
            vcf=call_j.output_vcf,
            shifted_vcf=shifted_call_j.output_vcf,
            reference=mito_ref,
            shift_back_chain=shifted_mito_ref.shift_back_chain,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_j)
        assert isinstance(merge_j.output_vcf, hb.ResourceGroup)

        # Merge the mutect stats output files (needed for filtering)
        merge_stats_J = mito.merge_mutect_stats(
            b=get_batch(),
            first_stats_file=call_j.output_vcf['vcf.gz.stats'],
            second_stats_file=shifted_call_j.output_vcf['vcf.gz.stats'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_stats_J)
        assert isinstance(merge_stats_J.combined_stats, hb.ResourceFile)

        # Initial round of filtering to exclude blacklist and high alt alleles
        initial_filter_j = mito.filter_variants(
            b=get_batch(),
            vcf=merge_j.output_vcf,
            reference=mito_ref,
            merged_mutect_stats=merge_stats_J.combined_stats,
            # alt_allele and vaf config hardcoded in this round of filtering as per
            # https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=4,
            min_allele_fraction=0,
            f_score_beta=get_config()['mito_snv']['f_score_beta'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(initial_filter_j)
        assert isinstance(initial_filter_j.output_vcf, hb.ResourceGroup)

        # SplitMultiAllelics AND remove non-passing sites
        # Output is only used for input to haplocheck
        split_multiallelics_j = mito.split_multi_allelics(
            b=get_batch(),
            vcf=initial_filter_j.output_vcf,
            reference=mito_ref,
            remove_non_pass_sites=True,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(split_multiallelics_j)
        assert isinstance(split_multiallelics_j.output_vcf, hb.ResourceGroup)

        # Estimate level of contamination from mito reads
        get_contamination_j = mito.get_contamination(
            b=get_batch(),
            vcf=split_multiallelics_j.output_vcf,
            haplocheck_output=self.expected_outputs(sequencing_group)[
                'haplocheck_metrics'
            ],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(get_contamination_j)
        assert isinstance(get_contamination_j.haplocheck_output, hb.ResourceFile)

        # Parse contamination estimate reports
        parse_contamination_j, contamination_level = mito.parse_contamination_results(
            b=get_batch(),
            haplocheck_output=get_contamination_j.haplocheck_output,
            verifybamid_output=verifybamid_output,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(parse_contamination_j)

        # Filter round 2 - remove variants with VAF below estimated contamination
        second_filter_j = mito.filter_variants(
            b=get_batch(),
            vcf=initial_filter_j.output_vcf,
            reference=mito_ref,
            merged_mutect_stats=merge_stats_J.combined_stats,
            # alt_allele config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=4,
            min_allele_fraction=get_config()['mito_snv']['vaf_filter_threshold'],
            f_score_beta=get_config()['mito_snv']['f_score_beta'],
            # contamination_estimate=get_contamination_j.max_contamination,
            contamination_estimate=contamination_level.as_str(),
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(second_filter_j)
        assert isinstance(second_filter_j.output_vcf, hb.ResourceGroup)

        # Generate final output vcf
        split_multiallelics_j = mito.split_multi_allelics(
            b=get_batch(),
            vcf=second_filter_j.output_vcf,
            reference=mito_ref,
            remove_non_pass_sites=False,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(split_multiallelics_j)

        # Write the final vcf to the bucket
        get_batch().write_output(
            split_multiallelics_j.output_vcf,
            str(self.expected_outputs(sequencing_group)['out_vcf']).replace(
                '.vcf.bgz', ''
            ),
        )

        return self.make_outputs(
            sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs
        )


@stage(required_stages=[RealignMito, GenotypeMito])
class JoinMito(CohortStage):
    """
    Join and annotate Mito snv and indel calls
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a matrix table.
        """
        return {
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.tmp_prefix),
            'coverage_ht': self.prefix / 'mito' / 'coverage.ht',
            'coverage_mt': self.prefix / 'mito' / 'coverage.mt',
            'cohort_mt': self.prefix / 'mito' / 'cohort.mt',
            'cohort_vcf': self.prefix / 'mito' / 'cohort.vcf.bgz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run hail query scripts
        """
        jobs = []

        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['tmp_prefix']) / 'checkpoints'
        )

        base_level_coverage_by_sgid = {
            sequencing_group.id: inputs.as_path(
                target=sequencing_group,
                stage=RealignMito,
                key='base_level_coverage_metrics',
            )
            for sequencing_group in cohort.get_sequencing_groups()
        }

        vcf_path_by_sgid = {
            sequencing_group.id: str(
                inputs.as_path(
                    target=sequencing_group, stage=GenotypeMito, key='out_vcf'
                )
            )
            for sequencing_group in cohort.get_sequencing_groups()
        }

        # Summerise per-sample coverage reports into cohort level statistics
        annotate_coverage_j = mito_cohort.annotate_coverage(
            b=get_batch(),
            base_level_coverage_by_sgid=base_level_coverage_by_sgid,
            coverage_ht=self.expected_outputs(cohort)['coverage_ht'],
            job_attrs=self.get_job_attrs(cohort),
            checkpoint_prefix=checkpoint_prefix,
        )
        jobs.append(annotate_coverage_j)

        # Merge individual vcfs into a single cohort mt
        combine_vcfs_j = mito_cohort.combine_vcfs(
            b=get_batch(),
            vcf_path_by_sgid=vcf_path_by_sgid,
            coverage_mt_path=self.expected_outputs(cohort)['coverage_mt'],
            artifact_prone_sites_path=to_path(
                'gs://cpg-common-test/mito/artifact_prone_sites.hg38.bed'
            ),
            combined_vcf_mt_path=self.expected_outputs(cohort)['cohort_mt'],
            checkpoint_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(cohort),
        )
        combine_vcfs_j.depends_on(annotate_coverage_j)
        jobs.append(combine_vcfs_j)

        return self.make_outputs(
            cohort,
            data=self.expected_outputs(cohort),
            jobs=jobs,
        )


@stage(required_stages=[RealignMito, GenotypeMito, JoinMito])
class VepMito(CohortStage):
    """
    Run VEP on Mito VCF
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a matrix table.
        """
        return {
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.tmp_prefix),
            'sites_only_vcf': self.prefix / 'mito' / 'cohort.sites_only.vcf.gz',
            'vep_json': self.prefix / 'mito' / 'cohort.sites_only.vep.json',
            'mito_vep_ht': self.prefix / 'mito' / 'cohort.sites_only.vep.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """ """
        jobs = []

        # Input resources
        cohort_vcf = get_batch().read_input_group(
            **{
                'vcf.gz': str(inputs.as_path(
                    target=cohort, stage=JoinMito, key='cohort_vcf'
                ))
            }
        )

        # Make a sites only vcf to run vep on
        siteonly_j, siteonly_vcf = joint_genotyping.add_make_sitesonly_job(
            b=get_batch(),
            input_vcf=cohort_vcf,
            output_vcf_path=self.expected_outputs(cohort)['sites_only_vcf'],
            job_attrs=self.get_job_attrs(cohort),
        )
        if siteonly_j:
            jobs.append(siteonly_j)

        # Run VEP
        vep_j = vep.vep_one(
            get_batch(),
            vcf=self.expected_outputs(cohort)['sites_only_vcf'],
            out_path=self.expected_outputs(cohort)['vep_json'],
            out_format='json',
            job_attrs=self.get_job_attrs(cohort),
        )
        if vep_j:
            if jobs:
                vep_j.depends_on(*jobs)
            jobs.append(vep_j)

        # Convert to VEP output to a hail table
        json_to_ht_j = vep.gather_vep_json_to_ht(
            get_batch(),
            vep_results_paths=[self.expected_outputs(cohort)['vep_json'],],
            out_path=self.expected_outputs(cohort)['mito_vep_ht'],
            job_attrs=self.get_job_attrs(cohort),
            depends_on=jobs
        )
        json_to_ht_j.depends_on(*jobs)
        jobs.append(json_to_ht_j)

        return self.make_outputs(
            cohort,
            data=self.expected_outputs(cohort),
            jobs=jobs,
        )


@stage(required_stages=[JoinMito, VepMito])
class AnotateMito(CohortStage):
    """
    Annotate final MT
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a matrix table.
        """
        return {
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.tmp_prefix),
            'mt': self.prefix / 'mito' / 'cohort.annotated.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """ """

        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['tmp_prefix']) / 'checkpoints'
        )

        # Make a sites only vcf to run vep on
        jobs = seqr_loader.annotate_cohort_jobs(
            b=get_batch(),
            vcf_path=inputs.as_path(cohort, JoinMito, 'cohort_mt'),
            out_mt_path=self.expected_outputs(cohort)['mt'],
            checkpoint_prefix=checkpoint_prefix,
            vep_ht_path=inputs.as_path(cohort, VepMito, 'mito_vep_ht'),
            job_attrs=self.get_job_attrs(cohort),
            use_dataproc=False
        )

        return self.make_outputs(
                    cohort,
                    data=self.expected_outputs(cohort),
                    jobs=jobs,
                )
