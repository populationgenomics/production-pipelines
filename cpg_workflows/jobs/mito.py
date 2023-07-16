"""
Hail Batch jobs needed to call mitochondrial SNVs
"""
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath


def subset_cram_to_chrM(
    b,
    cram_path: CramPath,
    job_attrs: dict | None = None,
) -> Job:
    """
    Extract read pairs with at least one read mapping to chrM

    Specific config for selection of which reads to extract may influence the impact of
    NUMTs. Because of this we stick to the exact command used by Broad until we have a
    good reason not to.

    Args:
        cram_path: Input cram to extract chrM reads from.

    Outputs:
        job.output_bam: a ResourceGroup containing bam of reads mapped to chrM

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

    # We are only accessing a tiny fraction of the genome. Mounting is the best option.
    bucket = cram_path.path.drive
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

    j.command(command(cmd))
    return j


def mito_realign(
    b,
    sequencing_group_id: str,
    input_bam: hb.ResourceGroup,
    mito_ref: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    Re-align reads to mitochondrial genome

    Uses bazam to extract fastq from input bam then pipes this directly to bwa for
    mapping.

    Args:
        sequencing_group_id: CPG sequencing_group id for inclusion in RG header
        input_bam: Bam for realignment
        mito_ref: Resource group containing the mito genome and bwa index to align to

    Outputs:
        job.output_cram: A sorted cram

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
            -R '@RG\\tID:{sequencing_group_id}\\tSM:{sequencing_group_id}' \
            - | \
        samtools view -bSu -T {mito_ref.base} - | \
        samtools sort -o {j.output_cram}
        """
    j.command(command(cmd))

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

    """

    j.command(command(cmd))
    if metrics:
        b.write_output(j.metrics, str(metrics))
    if metrics:
        b.write_output(j.theoretical_sensitivity, str(theoretical_sensitivity))

    return j


def extract_coverage_mean(
    b,
    metrics: hb.ResourceFile,
    mean_path: Path | None = None,
    median_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Extract mean and median coverage values for sequencing group

    Args:
        metrics: CollectWgsMetrics metrics output to process.
        mean_path: Output path for mean file.
        median_path: Output path for median file.

    Outputs:
        job.mean_coverage: mean coverage of chrM
        job.median_coverage: median coverage of chrM
    """
    job_attrs = job_attrs or {}
    j = b.new_job('extract_coverage_mean', job_attrs)
    j.image(image_path('peer'))


    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
        R --vanilla <<CODE
        df = read.table(
            "{metrics}",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\\\t',nrows=1
        )
        write.table(
            floor(df[,"MEAN_COVERAGE"]),
            "{j.mean_coverage}",
            quote=F, col.names=F, row.names=F)
        write.table(
            df[,"MEDIAN_COVERAGE"],
            "{j.median_coverage}",
            quote=F, col.names=F, row.names=F)
        CODE
    """

    j.command(command(cmd))
    if mean_path:
        b.write_output(j.mean_coverage, str(mean_path))
    if median_path:
        b.write_output(j.median_coverage, str(median_path))

    return j


def coverage_at_every_base(
    b,
    cram: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    intervals_list: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run picard CollectHsMetrics to calculate read coverage at each base.

    Args:
        cram: Input cram
        reference: Cram reference genome
        intervals_list: intervals list of target region

    Outputs:
        job.per_base_coverage
        job.hs_metics_out
    """
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
        job.merged_coverage: Merged coverage tsv.
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
      write.table(combined_table, "{j.merged_coverage}", row.names=F, col.names=T, quote=F, sep="\\\t")

    CODE
    """

    j.command(command(cmd))
    b.write_output(j.merged_coverage, str(merged_coverage))

    return j
