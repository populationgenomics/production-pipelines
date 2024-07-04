"""
Stage to call SNVs in the mitochondrial genome of a single sequencing_group.

Reimplemented version of;
https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl

"""

from functools import cache

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import mito, picard, vep
from cpg_workflows.stages.align import Align
from cpg_workflows.stages.cram_qc import CramQC
from cpg_workflows.workflow import (
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@cache
def get_mito_references(ref_path: str = 'gnomad_mito', shifted: bool = False) -> hb.ResourceGroup:
    """
    get various mito config entries, reads them into the current batch
    single method switches between shifted and non-shifted references
    Args:
        ref_path: path in config to the mito reference files
        shifted: whether to get the shifted reference files

    Returns:
        dict: mito config entries
    """
    shifted_str = 'shifted_' if shifted else ''
    mito_fa = reference_path(f'{ref_path}/{shifted_str}fasta')
    return get_batch().read_input_group(
        dict=reference_path(f'{ref_path}/{shifted_str}dict'),
        base=mito_fa,
        amb=mito_fa + '.amb',
        ann=mito_fa + '.ann',
        bwt=mito_fa + '.bwt',
        fai=mito_fa + '.fai',
        pac=mito_fa + '.pac',
        sa=mito_fa + '.sa',
    )


@cache
def get_control_region_intervals() -> hb.ResourceGroup:
    """
    get mito control region intervals
    Returns:
        dict: mito control region intervals
    """

    return get_batch().read_input_group(
        control_region_shifted=reference_path('gnomad_mito/shifted_control_region_interval'),
        non_control_region=reference_path('gnomad_mito/non_control_region_interval'),
    )


# alt_allele config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
MAX_ALT_ALLELE_COUNT = 4


@stage(
    required_stages=Align,
    analysis_type='mito-cram',
    analysis_keys=['non_shifted_cram'],
)
class RealignMito(SequencingGroupStage):
    """
    Re-align mitochondrial genome of a single sequencing_group.

    This is a re-implementation of the broad pipeline (as of 03/22) used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

    Mitochondrial variant calling can be subject to subtle artifacts resulting
    from mapping errors and other complexities such as NUMTs. In an attempt to minimize
    these issues, this pipeline attempts to re-implement as much of the logic, tools and
    configuration as possible.

    The main phases of analysis include:
        - Extraction of reads mapping to chrM from the main sequencing_group cram.
        - Realignment of chrM reads to chrM reference using bwa.
        - Realignment of chrM reads to a "shifted" chrM reference using bwa. The "shifted"
            reference is the same sequence but with a different linearisation point. This
            is used to overcome mapping artifacts at the boundary (aka the "control
            region" of the mt genome).
        - Generation of overall and  base level alignment statistics.

    Mitochondrial reference indexes and lift over chains were sourced from Broad
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
            'base_level_coverage_metrics': main / 'mito' / f'{sequencing_group.id}.base_level_coverage.tsv',
            'coverage_metrics': analysis / 'mito' / f'{sequencing_group.id}.coverage_metrics.txt',
            'coverage_mean': analysis / 'mito' / f'{sequencing_group.id}.coverage_mean.txt',
            'coverage_median': analysis / 'mito' / f'{sequencing_group.id}.coverage_median.txt',
            'theoretical_sensitivity_metrics': analysis / 'mito' / f'{sequencing_group.id}.theoretical_sensitivity.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # Mitochondrial specific reference files.
        mito_ref = get_mito_references()
        shifted_mito_ref = get_mito_references(shifted=True)
        intervals = get_control_region_intervals()

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
        assert isinstance(realign_mkdup_j, Job)
        jobs.append(realign_mkdup_j)
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
        assert isinstance(shifted_mkdup_j, Job)
        jobs.append(shifted_mkdup_j)
        assert isinstance(shifted_mkdup_j.output_cram, hb.ResourceGroup)

        # Collect coverage metrics (only on non-shifted)
        coverage_metrics_J = mito.collect_coverage_metrics(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            reference=mito_ref,
            metrics=self.expected_outputs(sequencing_group)['coverage_metrics'],
            theoretical_sensitivity=self.expected_outputs(sequencing_group)['theoretical_sensitivity_metrics'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(coverage_metrics_J)
        assert isinstance(coverage_metrics_J.metrics, hb.ResourceFile)

        # Extract mean and median coverage (only on non-shifted)
        extract_coverage_mean_j = mito.extract_coverage_mean(
            b=get_batch(),
            metrics=coverage_metrics_J.metrics,
            mean_path=self.expected_outputs(sequencing_group)['coverage_mean'],
            median_path=self.expected_outputs(sequencing_group)['coverage_median'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(extract_coverage_mean_j)

        # CoverageAtEveryBase
        # Provides coverage at each base so low coverage sites can be
        # considered ./. rather than 0/0.
        non_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=realign_mkdup_j.output_cram,
            intervals_list=intervals.non_control_region,
            reference=mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(non_control_region_coverage_j)
        assert isinstance(non_control_region_coverage_j.per_base_coverage, hb.ResourceFile)

        shifted_control_region_coverage_j = mito.coverage_at_every_base(
            b=get_batch(),
            cram=shifted_mkdup_j.output_cram,
            intervals_list=intervals.control_region_shifted,
            reference=shifted_mito_ref,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(shifted_control_region_coverage_j)
        assert isinstance(shifted_control_region_coverage_j.per_base_coverage, hb.ResourceFile)

        # Merge coverage stats
        merge_coverage_j = mito.merge_coverage(
            b=get_batch(),
            non_cr_coverage=non_control_region_coverage_j.per_base_coverage,
            shifted_cr_coverage=shifted_control_region_coverage_j.per_base_coverage,
            merged_coverage=self.expected_outputs(sequencing_group)['base_level_coverage_metrics'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_coverage_j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@stage(required_stages=[RealignMito, CramQC])
class GenotypeMito(SequencingGroupStage):
    """
    Call SNVs in the mitochondrial genome of a single sequencing_group.
    This is a re-implementation of the broad pipeline as of 03/22 that was used for
    gnomAD v3 and broad seqr:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl
    A default config file here:
    https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/mitochondria_m2_wdl/ExampleInputsMitochondriaPipeline.json

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
            'haplocheck_metrics': analysis / 'mito' / f'{sequencing_group.id}.haplocheck.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # Mitochondrial specific reference files.
        mito_ref = get_mito_references()
        shifted_mito_ref = get_mito_references(shifted=True)

        jobs = []

        # Get input resources
        non_shifted_cram = get_batch().read_input_group(
            cram=str(inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram')),
            crai=str(inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram')) + '.crai',
        )
        shifted_cram = get_batch().read_input_group(
            cram=str(inputs.as_path(sequencing_group, RealignMito, 'shifted_cram')),
            crai=str(inputs.as_path(sequencing_group, RealignMito, 'shifted_cram')) + '.crai',
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

        # read in the shift-back chain file
        shift_back_chain = get_batch().read_input(str(reference_path('gnomad_mito/shift_back_chain')))

        # Merge the wt and shifted VCFs
        merge_j = mito.liftover_and_combine_vcfs(
            b=get_batch(),
            vcf=call_j.output_vcf,
            shifted_vcf=shifted_call_j.output_vcf,
            reference=mito_ref,
            shift_back_chain=shift_back_chain,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_j)
        assert isinstance(merge_j.output_vcf, hb.ResourceGroup)

        # Merge the mutect stats output files (needed for filtering)
        merge_stats_j = mito.merge_mutect_stats(
            b=get_batch(),
            first_stats_file=call_j.output_vcf['vcf.gz.stats'],
            second_stats_file=shifted_call_j.output_vcf['vcf.gz.stats'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(merge_stats_j)
        assert isinstance(merge_stats_j.combined_stats, hb.ResourceFile)

        # Initial round of filtering to exclude blacklist and high alt alleles
        initial_filter_j = mito.filter_variants(
            b=get_batch(),
            vcf=merge_j.output_vcf,
            reference=mito_ref,
            merged_mutect_stats=merge_stats_j.combined_stats,
            # alt_allele and vaf config hardcoded in this round of filtering as per
            # https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=MAX_ALT_ALLELE_COUNT,
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
            haplocheck_output=self.expected_outputs(sequencing_group)['haplocheck_metrics'],
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
            merged_mutect_stats=merge_stats_j.combined_stats,
            # alt_allele config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
            max_alt_allele_count=MAX_ALT_ALLELE_COUNT,
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
        output_vcf_root = str(self.expected_outputs(sequencing_group)['out_vcf']).replace('.vcf.bgz', '')

        get_batch().write_output(
            split_multiallelics_j.output_vcf,
            output_vcf_root,
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)


@stage(
    required_stages=[RealignMito, GenotypeMito],
    analysis_type='web',
    analysis_keys=['mitoreport'],
)
class MitoReport(SequencingGroupStage):
    """
    Run the standalone MitoReport program on each SG

    This is not part of the Broad Mito pipeline, but generates an alternative non-seqr
    based interpretable html report of mito variants.

    Requires vep annotated individual vcf.

    https://github.com/bioinfomethods/mitoreport
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        main = sequencing_group.dataset.prefix()
        web = sequencing_group.dataset.web_prefix()
        return {
            'vep_vcf': main / 'mito' / f'{sequencing_group.id}.mito.vep.vcf.gz',
            'mitoreport': web / 'mito' / f'mitoreport-{sequencing_group.id}' / 'index.html',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        mito_ref = get_mito_references()
        jobs = []

        vep_j = vep.vep_one(
            get_batch(),
            vcf=inputs.as_path(sequencing_group, GenotypeMito, 'out_vcf'),
            out_path=self.expected_outputs(sequencing_group)['vep_vcf'],
            out_format='vcf',
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        if vep_j:
            jobs.append(vep_j)

        mitoreport_j = mito.mitoreport(
            get_batch(),
            sequencing_group=sequencing_group,
            vcf_path=self.expected_outputs(sequencing_group)['vep_vcf'],
            cram_path=inputs.as_path(sequencing_group, RealignMito, 'non_shifted_cram'),
            mito_ref=mito_ref,
            output_path=self.expected_outputs(sequencing_group)['mitoreport'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        if mitoreport_j:
            mitoreport_j.depends_on(*jobs)
            jobs.append(mitoreport_j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
