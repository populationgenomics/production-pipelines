"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from cpg_utils.hail_batch import command
from cpg_workflows.resources import HIGHMEM, STANDARD, storage_for_cram_qc_job
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse, exists


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    source_intervals_path: Path | None = None,
    job_attrs: dict[str, str] | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job | None, list[hb.ResourceFile]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param source_intervals_path: path to source intervals to split. Would check for
        config if not provided.
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredictable number of intervals. WDL can
    handle that, but Hail Batch is not dynamic and expects a certain number
    of output files.
    """
    assert scatter_count > 0, scatter_count
    sequencing_type = get_config()['workflow']['sequencing_type']
    source_intervals_path = source_intervals_path or reference_path(
        f'broad/{sequencing_type}_calling_interval_lists'
    )

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(str(source_intervals_path))]

    if output_prefix and exists(output_prefix / '1.interval_list'):
        return None, [
            b.read_input(str(output_prefix / f'{idx + 1}.interval_list'))
            for idx in range(scatter_count)
        ]

    if not source_intervals_path and exists(
        (
            existing_split_intervals_prefix := (
                reference_path('intervals_prefix')
                / sequencing_type
                / f'{scatter_count}intervals'
            )
        )
        / '1.interval_list'
    ):
        # We already have split intervals for this sequencing_type:
        return None, [
            b.read_input(
                str(existing_split_intervals_prefix / f'{idx + 1}.interval_list')
            )
            for idx in range(scatter_count)
        ]

    j = b.new_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | dict(tool='picard IntervalListTools'),
    )
    j.image(image_path('picard'))
    STANDARD.set_resources(j, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        'genome': 100000,
        'exome': 0,
    }.get(sequencing_type, 0)

    cmd = f"""
    mkdir $BATCH_TMPDIR/out
    
    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(source_intervals_path))} \
    OUTPUT=$BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln $BATCH_TMPDIR/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )

    intervals: list[hb.ResourceFile] = []
    for idx in range(scatter_count):
        interval = j[f'{idx + 1}.interval_list']
        assert isinstance(interval, hb.ResourceFile)
        intervals.append(interval)
    return j, intervals


def markdup(
    b: hb.Batch,
    sorted_bam: hb.ResourceFile,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    out_markdup_metrics_path: Path | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Make job that runs Picard MarkDuplicates and converts the result to CRAM.
    """
    job_attrs = (job_attrs or {}) | dict(tool='picard_MarkDuplicates')
    j = b.new_job('MarkDuplicates', job_attrs)
    if can_reuse(output_path, overwrite):
        j.name = f'{j.name} [reuse]'
        return j

    j.image(image_path('picard_samtools'))
    resource = HIGHMEM.request_resources(ncpu=4)
    # enough for input BAM and output CRAM
    resource.attach_disk_storage_gb = 250
    resource.set_to_job(j)
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )
    fasta_reference = fasta_res_group(b)

    assert isinstance(j.output_cram, hb.ResourceGroup)
    cmd = f"""
    picard MarkDuplicates -Xms{resource.get_java_mem_mb()}M \\
    I={sorted_bam} O=/dev/stdout M={j.markdup_metrics} \\
    TMP_DIR=$(dirname {j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate \\
    | samtools view -@{resource.get_nthreads() - 1} -T {fasta_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@{resource.get_nthreads() - 1} {j.output_cram.cram} {j.output_cram['cram.crai']}
    """
    j.command(command(cmd, monitor_space=True))

    if output_path:
        b.write_output(j.output_cram, str(output_path.with_suffix('')))
        if out_markdup_metrics_path:
            b.write_output(
                j.markdup_metrics,
                str(out_markdup_metrics_path),
            )

    return j


def vcf_qc(
    b: hb.Batch,
    vcf_or_gvcf: hb.ResourceGroup,
    is_gvcf: bool,
    job_attrs: dict | None = None,
    output_summary_path: Path | None = None,
    output_detail_path: Path | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Make job that runs Picard CollectVariantCallingMetrics.
    """
    if output_summary_path and can_reuse(output_detail_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard CollectVariantCallingMetrics'}
    j = b.new_job('CollectVariantCallingMetrics', job_attrs)
    j.image(image_path('picard'))
    res = STANDARD.set_resources(j, storage_gb=20, mem_gb=3)
    reference = fasta_res_group(b)
    dbsnp_vcf = b.read_input_group(
        base=str(reference_path('broad/dbsnp_vcf')),
        index=str(reference_path('broad/dbsnp_vcf_index')),
    )
    sequencing_type = get_config()['workflow']['sequencing_type']
    intervals_file = b.read_input(
        str(reference_path(f'broad/{sequencing_type}_evaluation_interval_lists'))
    )

    if is_gvcf:
        input_file = vcf_or_gvcf['g.vcf.gz']
    else:
        input_file = vcf_or_gvcf['vcf.gz']

    cmd = f"""\
    picard -Xms2000m -Xmx{res.get_java_mem_mb()}m \
    CollectVariantCallingMetrics \
    INPUT={input_file} \
    OUTPUT=$BATCH_TMPDIR/prefix \
    DBSNP={dbsnp_vcf['base']} \
    SEQUENCE_DICTIONARY={reference['dict']} \
    TARGET_INTERVALS={intervals_file} \
    GVCF_INPUT={"true" if is_gvcf else "false"}
    
    cp $BATCH_TMPDIR/prefix.variant_calling_summary_metrics {j.summary}
    cp $BATCH_TMPDIR/prefix.variant_calling_detail_metrics {j.detail}
    """

    j.command(command(cmd))

    if output_summary_path:
        b.write_output(j.summary, str(output_summary_path))
    if output_detail_path:
        b.write_output(j.detail, str(output_detail_path))
    return j


def picard_collect_metrics(
    b,
    cram_path: CramPath,
    out_alignment_summary_metrics_path: Path,
    out_base_distribution_by_cycle_metrics_path: Path,
    out_insert_size_metrics_path: Path,
    out_quality_by_cycle_metrics_path: Path,
    out_quality_yield_metrics_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run picard CollectMultipleMetrics metrics for sample QC.
    Based on https://github.com/broadinstitute/warp/blob/master/tasks/broad/Qc.wdl#L141
    """
    if can_reuse(
        [
            out_alignment_summary_metrics_path,
            out_base_distribution_by_cycle_metrics_path,
            out_insert_size_metrics_path,
            out_quality_by_cycle_metrics_path,
            out_quality_yield_metrics_path,
        ],
        overwrite,
    ):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard_CollectMultipleMetrics'}
    j = b.new_job('Picard CollectMultipleMetrics', job_attrs)
    j.image(image_path('picard'))
    res = STANDARD.request_resources(ncpu=2)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    reference = fasta_res_group(b)

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    picard -Xmx{res.get_java_mem_mb()}m \\
      CollectMultipleMetrics \\
      INPUT=$CRAM \\
      REFERENCE_SEQUENCE={reference.base} \\
      OUTPUT=$BATCH_TMPDIR/prefix \\
      ASSUME_SORTED=true \\
      PROGRAM=null \\
      VALIDATION_STRINGENCY=SILENT \\
      PROGRAM=CollectAlignmentSummaryMetrics \\
      PROGRAM=CollectInsertSizeMetrics \\
      PROGRAM=MeanQualityByCycle \\
      PROGRAM=CollectBaseDistributionByCycle \\
      PROGRAM=CollectQualityYieldMetrics \\
      METRIC_ACCUMULATION_LEVEL=null \\
      METRIC_ACCUMULATION_LEVEL=SAMPLE
      
    ls $BATCH_TMPDIR/
    cp $BATCH_TMPDIR/prefix.alignment_summary_metrics {j.out_alignment_summary_metrics}
    cp $BATCH_TMPDIR/prefix.base_distribution_by_cycle_metrics {j.out_base_distribution_by_cycle_metrics}
    cp $BATCH_TMPDIR/prefix.insert_size_metrics {j.out_insert_size_metrics}
    cp $BATCH_TMPDIR/prefix.quality_by_cycle_metrics {j.out_quality_by_cycle_metrics}
    cp $BATCH_TMPDIR/prefix.quality_yield_metrics {j.out_quality_yield_metrics}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(
        j.out_alignment_summary_metrics, str(out_alignment_summary_metrics_path)
    )
    b.write_output(j.out_insert_size_metrics, str(out_insert_size_metrics_path))
    b.write_output(
        j.out_quality_by_cycle_metrics, str(out_quality_by_cycle_metrics_path)
    )
    b.write_output(
        j.out_base_distribution_by_cycle_metrics,
        str(out_base_distribution_by_cycle_metrics_path),
    )
    b.write_output(j.out_quality_yield_metrics, str(out_quality_yield_metrics_path))
    return j


def picard_hs_metrics(
    b,
    cram_path: CramPath,
    job_attrs: dict | None = None,
    out_picard_hs_metrics_path: Path | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run picard CollectHsMetrics metrics.
    Based on https://github.com/broadinstitute/warp/blob/master/tasks/broad/Qc.wdl#L528
    """
    if can_reuse(out_picard_hs_metrics_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard_CollectHsMetrics'}
    j = b.new_job('Picard CollectHsMetrics', job_attrs)
    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'exome'
    res = STANDARD.request_resources(ncpu=2)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    reference = fasta_res_group(b)
    interval_file = b.read_input(
        str(reference_path('broad/exome_evaluation_interval_lists'))
    )

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    # Picard is strict about the interval-list file header - contigs md5s, etc. - and 
    # if md5s do not match the ref.dict file, picard would crash. So fixing the header 
    # by converting the interval-list to bed (i.e. effectively dropping the header) 
    # and back to interval-list (effectively re-adding the header from input ref-dict).
    # VALIDATION_STRINGENCY=SILENT does not help.
    picard IntervalListToBed \\
    I={interval_file} \\
    O=$BATCH_TMPDIR/intervals.bed
    picard BedToIntervalList \\
    I=$BATCH_TMPDIR/intervals.bed \\
    O=$BATCH_TMPDIR/intervals.interval_list \\
    SD={reference.dict}
    
    picard -Xmx{res.get_java_mem_mb()}m \\
      CollectHsMetrics \\
      INPUT=$CRAM \\
      REFERENCE_SEQUENCE={reference.base} \\
      VALIDATION_STRINGENCY=SILENT \\
      TARGET_INTERVALS=$BATCH_TMPDIR/intervals.interval_list \\
      BAIT_INTERVALS=$BATCH_TMPDIR/intervals.interval_list \\
      METRIC_ACCUMULATION_LEVEL=null \\
      METRIC_ACCUMULATION_LEVEL=SAMPLE \\
      METRIC_ACCUMULATION_LEVEL=LIBRARY \\
      OUTPUT={j.out_hs_metrics}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.out_hs_metrics, str(out_picard_hs_metrics_path))
    return j


def picard_wgs_metrics(
    b,
    cram_path: CramPath,
    out_picard_wgs_metrics_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
    read_length: int = 250,
) -> Job | None:
    """
    Run picard CollectWgsMetrics metrics.
    Based on https://github.com/broadinstitute/warp/blob/e1ac6718efd7475ca373b7988f81e54efab608b4/tasks/broad/Qc.wdl#L444
    """
    if can_reuse(out_picard_wgs_metrics_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard_CollectWgsMetrics'}
    j = b.new_job('Picard CollectWgsMetrics', job_attrs)

    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'genome'
    res = STANDARD.request_resources(ncpu=2)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    reference = fasta_res_group(b)
    interval_file = b.read_input(
        str(reference_path('broad/genome_coverage_interval_list'))
    )

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    picard -Xmx{res.get_java_mem_mb()}m \\
      CollectWgsMetrics \\
      INPUT=$CRAM \\
      VALIDATION_STRINGENCY=SILENT \\
      REFERENCE_SEQUENCE={reference.base} \\
      INTERVALS={interval_file} \\
      OUTPUT={j.out_csv} \\
      USE_FAST_ALGORITHM=true \\
      READ_LENGTH={read_length}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.out_csv, str(out_picard_wgs_metrics_path))
    return j
