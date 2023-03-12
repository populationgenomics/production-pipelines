"""
Create Hail Batch jobs to call mitochondrial SNVs
"""
import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.resource import PythonResult

from cpg_utils.config import get_config
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse

from cpg_workflows.jobs import picard

# from cpg_workflows.mito_pipeline_scripts import get_final_contamination


def subset_cram_to_chrM(
    b,
    cram_path: CramPath,
    job_attrs: dict | None = None,
) -> Job:
    """
    Extract read pairs with at least one read mapping to chrM

    Specific config for selection of which reads to extract may influence the impact of
    NUMTs. Because of this we stick to the exact command used by Broad until we know
    better.

    Args:
        cram_path: Input cram to extract chrM reads from.

    Outputs:
        output_bam: a ResourceGroup containing bam of reads mapped to chrM

    Cmd from:
    https://github.com/broadinstitute/gatk/blob/330c59a5bcda6a837a545afd2d453361f373fae3/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl#LL188-L244C2

    """

    job_attrs = job_attrs or {}
    j = b.new_job('subset_cram_to_chrM', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    reference = fasta_res_group(b)
    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bai',
        }
    )

    # We are only accessing a fraction fof the genome. Mounting is the best option.
    bucket = cram_path.path.drive
    print(f'bucket = {bucket}')
    bucket_mount_path = to_path('/bucket')
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)
    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])

    cmd = f"""

        gatk PrintReads \
            -R {reference.base} \
            -L chrM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -I {mounted_cram_path} \
            -O {j.output_bam.bam}

    """

    j.command(command(cmd, define_retry_function=True))
    return j


def mito_realign(
    b,
    sample_id: str,
    input_bam: hb.ResourceGroup,
    mito_ref: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    Re-align reads to mitochondrial genome

    Uses bazam to extract fastq from input bam then pipes this directly to bwa for
    mapping.

    Args:
        sample_id: CPG sample id for inclusion in RG header
        input_bam: Bam for realignment

    Outputs:
        output_cram: A sorted cram

    Bwa command from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignmentPipeline.wdl#L59

    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_realign', job_attrs)
    j.image(image_path('bwa'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    nthreads = res.get_nthreads()

    cmd = f"""\
        bazam -Xmx16g \
            -n{min(nthreads, 6)} -bam {input_bam.bam} | \
        bwa \
            mem -K 100000000 -p -v 3 -t 2 -Y {mito_ref.base} \
            -R '@RG\\tID:{sample_id}\\tSM:{sample_id}' \
            - | \
        samtools view -bSu -T {mito_ref.base} - | \
        samtools sort -o {j.output_cram}
        """
    j.command(command(cmd, define_retry_function=True))

    return j


def collect_coverage_metrics(
    b,
    cram: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    metrics: Path | None = None,
    theoretical_sensitivity: Path | None = None,
    read_length_for_optimization: int = 151,
    coverage_cap: int = 100000,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run CollectWgsMetrics

    Args:
        cram: Input Cram
        reference: reference fastq
        read_length_for_optimization:  Read length used for optimization only. If this is
            too small CollectWgsMetrics might fail, but the results are not affected
            by this number. [Default: 151] [default: 100000 (from wdl)].
        coverage_cap: Treat positions with coverage exceeding this value as if they had
            coverage at this value.

    Outputs:
        metrics: output file
        theoretical_sensitivity: Undocumented CollectWgsMetrics output?
    """
    job_attrs = job_attrs or {}
    j = b.new_job('collect_coverage_metrics', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
        picard \
            CollectWgsMetrics \
            INPUT={cram.cram} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE={reference.base} \
            OUTPUT={j.metrics} \
            USE_FAST_ALGORITHM=true \
            READ_LENGTH={read_length_for_optimization} \
            COVERAGE_CAP={coverage_cap} \
            INCLUDE_BQ_HISTOGRAM=true \
            THEORETICAL_SENSITIVITY_OUTPUT={j.theoretical_sensitivity}

            # TODO: build an picard image with R to extract values
            # R --vanilla <<CODE
            # df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
            # write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
            # write.table(df[,"MEDIAN_COVERAGE"], "median_coverage.txt", quote=F, col.names=F, row.names=F)
            # CODE
    """

    j.command(command(cmd))
    if metrics:
        b.write_output(j.metrics, str(metrics))
    if metrics:
        b.write_output(j.theoretical_sensitivity, str(theoretical_sensitivity))

    return j


def coverage_at_every_base(
    b,
    cram: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    intervals_list: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """ """
    job_attrs = job_attrs or {}
    j = b.new_job('coverage_at_every_base', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
    picard CollectHsMetrics \
      I={cram.cram} \
      R={reference.base} \
      PER_BASE_COVERAGE={j.per_base_coverage} \
      O={j.hs_metics_out} \
      TI={intervals_list} \
      BI={intervals_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    """

    j.command(command(cmd))

    return j


def merge_coverage(
    b,
    non_cr_coverage: hb.ResourceFile,
    shifted_cr_coverage: hb.ResourceFile,
    merged_coverage: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Merge per base coverage files from non-control region and shifted control region.

    Args:
        non_cr_coverage: Per-base coverage from CollectHsMetrics for all of
            chrM excluding the control region.
        shifted_cr_coverage: Per-base coverage from CollectHsMetrics for only the chrM
            control region in shifted coord space.
        merged_coverage: Path to write merged coverage tsv.

    Outputs:
        merged_coverage: Merged coverage tsv.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('merge_coverage', job_attrs)
    j.image(image_path('peer'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
    R --vanilla <<CODE
      shift_back = function(x) {{
        if (x < 8570) {{
          return(x + 8000)
        }} else {{
          return (x - 8569)
        }}
      }}

      control_region_shifted = read.table("{shifted_cr_coverage}", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("{non_cr_coverage}", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      write.table(combined_table, "{j.merged_coverage}", row.names=F, col.names=T, quote=F, sep="\t")

    CODE
    """

    j.command(command(cmd))
    b.write_output(j.merged_coverage, str(merged_coverage))

    return j


def mito_mutect2(
    b,
    cram: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    region: str,
    max_reads_per_alignment_start: int = 75,
    job_attrs: dict | None = None,
) -> Job:
    """
    Call SNPs and indels in mitochondrial genome using Mutect2 in"mitochondria-mode"
    Args:
        cram: Cam to call variants in.
        reference: Resource group of reference sequence to align to.
        region: Coordinate string restricting the region to call variants within.
        max_reads_per_alignment_start: Mutect argument. [Default: 75].
    Output:
        output_vcf: resource file containing vcf, index AND statistics file.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L417-L484
    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_mutect2', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    java_mem_mb = res.get_java_mem_mb()

    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            'vcf.gz.stats': '{root}.vcf.gz.stats',
        }
    )

    cmd = f"""
        gatk --java-options "-Xmx{java_mem_mb}m" Mutect2 \
            -R {reference.base} \
            -I {cram.cram} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {j.output_vcf['vcf.gz']} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            -L {region}
    """

    j.command(command(cmd))

    return j


def liftover_and_combine_vcfs(
    b,
    vcf: hb.ResourceGroup,
    shifted_vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    shift_back_chain: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """
    Lifts over shifted vcf and combines it with the rest of the chrM calls.
    Args:
        vcf: chrM variants mapped to wt chrM (exl control region).
        shifted_vcf: chrM control region variants mapped to shifted chrM.
        reference: resource group for wt chrM reference.
        shift_back_chain: chain file provided with shifted genome.
    Outputs:
        output_vcf: Merged vcf.gz in standard hg38 coordinate space.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL360-L415C2
    """
    job_attrs = job_attrs or {}
    j = b.new_job('liftover_and_combine_vcfs', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        lifted_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        picard LiftoverVcf \
            I={shifted_vcf['vcf.gz']} \
            O={j.lifted_vcf['vcf.gz']} \
            R={reference.base} \
            CHAIN={shift_back_chain} \
            REJECT={j.rejected_vcf}.vcf.gz
        picard MergeVcfs \
            I={vcf['vcf.gz']} \
            I={j.lifted_vcf['vcf.gz']} \
            O={j.output_vcf['vcf.gz']}
    """

    j.command(command(cmd))

    return j


def merge_mutect_stats(
    b,
    first_stats_file: hb.ResourceFile,
    second_stats_file: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """
    Merge stats files from two mutect runs
    Args:
        first_stats_file: Mutect stats ResourceFile
        second_stats_file: Mutect stats ResourceFile
    Outputs:
        combined_stats: Combined stats file
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL573-L598C2
    """
    job_attrs = job_attrs or {}
    j = b.new_job('merge_stats', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        gatk MergeMutectStats \
            --stats {first_stats_file} \
            --stats {second_stats_file} -O {j.combined_stats}
    """

    j.command(command(cmd))

    return j


def filter_variants(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    merged_mutect_stats: hb.ResourceFile,
    max_alt_allele_count: int,
    min_allele_fraction: int,
    contamination_estimate: hb.ResourceFile | None = None,
    f_score_beta: float = 1.0,
    job_attrs: dict | None = None,
) -> Job:
    """
    Filter mutect variant calls
    Uses:
      gatk FilterMutectCalls to filter on variant properties
      gatk VariantFiltration to exclude a (broad supplied) black list of problem variants.
    Args:
        vcf: input vcf
        merged_mutect_stats: Mutect statistics file
        max_alt_allele_count: Maximum alt alleles per site (VariantFiltration).
        min_allele_fraction: Hard cutoff for minimum allele fraction. All sites with VAF
            less than this cutoff will be filtered. (VariantFiltration).
        contamination_estimate: Estimated sample contamination level (VariantFiltration).
            This is passed as a single float in a file as it is calculated at an earlier
            step in the pipeline.
        f_score_beta: F score beta, the relative weight of recall to precision,
            used if OPTIMAL_F_SCORE strategy is chosen (VariantFiltration). Default: 1.0
            is default provided by VariantFiltration.
    Outputs:
        output_vcf: filtered vcf file.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL486-L571C2
    Note:
        contamination_estimate pre-calculation has been moved out of this function.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('filter_variants', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    blacklisted_sites = b.read_input_group(
        bed='gs://cpg-common-main/references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed',
        idx='gs://cpg-common-main/references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed.idx',
    )

    j.declare_resource_group(
        filtered_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    if contamination_estimate:
        contamination_estimate_string = (
            f'--contamination-estimate  $(cat {contamination_estimate})'
        )
    else:
        contamination_estimate_string = ''

    cmd = f"""
      gatk --java-options "-Xmx2500m" FilterMutectCalls -V {vcf['vcf.gz']} \\
        -R {reference.base} \\
        -O {j.filtered_vcf['vcf.gz']} \\
        --stats {merged_mutect_stats} \\
        --max-alt-allele-count {max_alt_allele_count} \\
        --mitochondria-mode \\
        --min-allele-fraction  {min_allele_fraction} \\
        --f-score-beta  {f_score_beta} \\
        {contamination_estimate_string}
      gatk VariantFiltration -V {j.filtered_vcf['vcf.gz']} \\
        -O {j.output_vcf['vcf.gz']} \\
        --apply-allele-specific-filters \\
        --mask {blacklisted_sites.bed} \\
        --mask-name "blacklisted_site"
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def split_multi_allelics(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    remove_non_pass_sites: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Splits multi allelics and removes non pass sites
    Uses LeftAlignAndTrimVariants to split then optionally use SelectVariants to select
    only passing variants.
    Args:
        vcf: Input vcf file.
        reference: chrM reference fasta.
        remove_non_pass_sites:
    Output:
        output_vcf: Final vcf file.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L600
    """
    job_attrs = job_attrs or {}
    j = b.new_job('split_multi_allelics', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(split_vcf={'vcf.gz': '{root}.vcf.gz'})
    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz'})

    cmd = f"""
        gatk LeftAlignAndTrimVariants \
            -R {reference.base} \
            -V {vcf['vcf.gz']} \
            -O {j.split_vcf['vcf.gz']} \
            --split-multi-allelics \
            --dont-trim-alleles \
            --keep-original-ac
        """
    if remove_non_pass_sites:
        cmd += f"""
            gatk SelectVariants \
                -V {j.split_vcf['vcf.gz']} \
                -O {j.output_vcf['vcf.gz']} \
                --exclude-filtered
        """
    else:
        cmd += f"""
            mv {j.split_vcf['vcf.gz']} {j.output_vcf['vcf.gz']}
        """

    j.command(command(cmd, define_retry_function=True))

    return j


def get_contamination(
    b,
    vcf: hb.ResourceGroup,
    haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Uses new HaplocheckerCLI to estimate levels of contamination in mitochondria based
    on known mitochondrial haplotypes.
    Args:
        vcf: input vcf of passing variants with multi-allelics split.
        haplocheck_output: Path to write results file.
    Outputs:
        max_contamination: ResourceFile containing a single float of final contamination
            estimate for sample.
    Cmd from:
      https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L239
    """
    job_attrs = job_attrs or {}
    j = b.new_job('get_contamination', job_attrs)
    j.image(image_path('haplocheckcli'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
        PARENT_DIR="$(dirname "{vcf['vcf.gz']}")"
        java -jar /haplocheckCLI.jar "$PARENT_DIR"
        cp output {j.haplocheck_output}

        # sed 's/\"//g' output > output-noquotes
        # grep "SampleID" output-noquotes > headers
        # FORMAT_ERROR="Bad contamination file format"
        # if [ `awk '{{print $2}}' headers` != "Contamination" ]; then
        #     echo $FORMAT_ERROR; exit 1
        # fi
        # if [ `awk '{{print $6}}' headers` != "HgMajor" ]; then
        #     echo $FORMAT_ERROR; exit 1
        # fi
        # if [ `awk '{{print $8}}' headers` != "HgMinor" ]; then
        #     echo $FORMAT_ERROR; exit 1
        # fi
        # if [ `awk '{{print $14}}' headers` != "MeanHetLevelMajor" ]; then
        #     echo $FORMAT_ERROR; exit 1
        # fi
        # if [ `awk '{{print $15}}' headers` != "MeanHetLevelMinor" ]; then
        #     echo $FORMAT_ERROR; exit 1
        # fi
        # # Extract values of interest from report (This feels nuts)
        # grep -v "SampleID" output-noquotes > output-data
        # has_contamination=$(awk -F "\t" '{{print $2}}' output-data)
        # major_level=$(awk -F "\t" '{{print $14}}' output-data)
        # minor_level=$(awk -F "\t" '{{print $15}}' output-data)
        # # Boil them down into a single value to use in filter
        # # Float arithmetic using strings... what could go wrong?
        # if [ $has_contamination == "YES" ]
        # then
        #     if [ $major_level == "0.0" ]
        #     then
        #         max_contamination=$minor_level
        #     else
        #         max_contamination=$(echo "1.0 - $major_level" | bc)
        #     fi
        # else
        #     max_contamination=0.0
        # fi
        # # Save to a file to pass to filter
        # echo "Estimated contamination = $max_contamination"
        # echo $max_contamination > {j.max_contamination}
        """

    j.command(command(cmd, define_retry_function=True))
    if haplocheck_output:
        b.write_output(j.haplocheck_output, str(haplocheck_output))

    return j


def parse_contamination_results(
    b,
    haplocheck_output: hb.ResourceFile,
    verifybamid_output: hb.ResourceFile | None = None,
    job_attrs: dict | None = None,
) -> tuple[Job, PythonResult]:
    """
    Try (again) to get haplocheck_output post processing done with a python job...

    Inputs:
        haplocheck_report: native output from haplocheckCLI
        verifybamid_output: [optional] native output from verifyBamID

        Based on logic here:
        https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#LL523-L524


    """
    job_attrs = job_attrs or {}
    j = b.new_python_job('parse_contamination_results', job_attrs)
    j.image(get_config()['workflow']['driver_image'])

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    # Hmmmm. So, the only way I can get this function to the execution node seems
    # to be by defining it locally within the job function. If I move it out to the same
    # scope as the job funtion (or anywhere else) hail seems to try to find it in the
    #  version of production_pipelines installed in the node (and it is not there unless
    # I publish a new image).
    # Someone please tell me there is a better way?
    def parse_contamination_worker(
        haplocheck_report: str, verifybamid_report: str | None
    ) -> float:
        """
        Process haplocheckCLI and verifyBamID outputs to get contamination level as a
        single float.
        """
        cleaned_lines = []
        with open(haplocheck_report) as haplocheck:
            # Split line on tabs and strip double quotes
            for line in haplocheck:
                cleaned_lines.append([x.strip('"') for x in line.strip().split('\t')])
        # sanity check and reformat
        assert len(cleaned_lines) == 2, "haplocheck report is unexpected format"
        assert len(cleaned_lines[0]) == 17, "haplocheck report is unexpected format"
        report = dict(zip(cleaned_lines[0], cleaned_lines[1]))

        # Determine final contamination level
        if report['Contamination'] == 'YES':
            if float(report['MeanHetLevelMajor']) == 0:
                max_contamination = float(report['MeanHetLevelMinor'])
            else:
                max_contamination = 1.0 - float(report['MeanHetLevelMajor'])
        else:
            max_contamination = 0.0

        # If verifybamid_report is provided, chose the higher of the two
        if verifybamid_report:
            with open(verifybamid_report) as verifybamid:
                lines = [line.split('\t') for line in verifybamid.readlines()]
                assert len(lines) == 2
                report = dict(zip(lines[0], lines[1]))

            verifybamid_estimate = float(report['FREEMIX'])

            if verifybamid_estimate > max_contamination:
                max_contamination = verifybamid_estimate

        return max_contamination

    contamination_level = j.call(
        parse_contamination_worker, haplocheck_output, verifybamid_output
    )

    return j, contamination_level
