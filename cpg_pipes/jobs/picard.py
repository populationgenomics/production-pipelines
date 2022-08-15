"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import HIGHMEM, STANDARD, storage_for_cram_qc_job
from cpg_pipes.filetypes import CramPath
from cpg_pipes.utils import can_reuse, exists


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    intervals_path: Path | str | None = None,
    job_attrs: dict | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job | None, list[hb.Resource | None]]:
    """
    Add a job that split genome into intervals to parallelise variant calling.

    As input interval file, takes intervals_path if provided, otherwise checks refs
    for the intervals of provided sequencing_type.

    This job calls picard's IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredicted number of intervals. WDL can handle
    that, but Hail Batch is not dynamic and have to expect certain number of output
    files.
    """
    assert scatter_count > 0

    sequencing_type = get_config()['workflow']['sequencing_type']

    if not intervals_path:
        # Did we cache split intervals for this sequencing_type?
        cache_bucket = (
            reference_path('intervals_prefix')
            / sequencing_type
            / f'{scatter_count}intervals'
        )
        if exists(cache_bucket / '1.interval_list'):
            return None, [
                b.read_input(str(cache_bucket / f'{idx + 1}.interval_list'))
                for idx in range(scatter_count)
            ]
        # Taking intervals file for the sequencing_type.
        intervals_path = reference_path(
            f'broad/{sequencing_type}_calling_interval_lists',
        )

    if scatter_count == 1:
        return None, [b.read_input(str(intervals_path))]

    job_attrs = (job_attrs or {}) | dict(tool='picard_IntervalListTools')
    job_name = f'Make {scatter_count} intervals for {sequencing_type}'

    if output_prefix and exists(output_prefix / '1.interval_list'):
        job_attrs['reuse'] = True
        return b.new_job(job_name, job_attrs), [
            b.read_input(str(output_prefix / f'{idx + 1}.interval_list'))
            for idx in range(scatter_count)
        ]

    j = b.new_job(job_name, job_attrs)
    j.image(image_path('picard'))
    STANDARD.set_resources(j, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        'genome': 100000,
        'exome': 0,
    }.get(sequencing_type, 0)

    cmd = f"""
    mkdir /io/batch/out
    
    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(intervals_path))} \
    OUTPUT=/io/batch/out
    ls /io/batch/out
    ls /io/batch/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln /io/batch/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(wrap_command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )
    return j, [j[f'{idx + 1}.interval_list'] for idx in range(scatter_count)]


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

    cmd = f"""
    picard MarkDuplicates -Xms{resource.get_java_mem_mb()}M \\
    I={sorted_bam} O=/dev/stdout M={j.markdup_metrics} \\
    TMP_DIR=$(dirname {j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate \\
    | samtools view -@{resource.get_nthreads() - 1} -T {fasta_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@{resource.get_nthreads() - 1} {j.output_cram.cram} {j.output_cram['cram.crai']}
    """
    j.command(wrap_command(cmd, monitor_space=True))

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

    job_attrs = (job_attrs or {}) | dict(tool='picard_CollectVariantCallingMetrics')
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
    OUTPUT=/io/batch/prefix \
    DBSNP={dbsnp_vcf['base']} \
    SEQUENCE_DICTIONARY={reference['dict']} \
    TARGET_INTERVALS={intervals_file} \
    GVCF_INPUT={"true" if is_gvcf else "false"}
    
    cp /io/batch/prefix.variant_calling_summary_metrics {j.summary}
    cp /io/batch/prefix.variant_calling_detail_metrics {j.detail}
    """

    j.command(wrap_command(cmd))

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

    job_attrs = (job_attrs or {}) | {'tool': 'picard'}
    j = b.new_job('Picard CollectMultipleMetrics', job_attrs)
    j.image(image_path('picard'))
    res = STANDARD.request_resources(mem_gb=7)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)

    cmd = f"""\
    picard -Xms5g -Xmx{res.get_java_mem_mb()}m \
      CollectMultipleMetrics \
      INPUT={cram.cram} \
      REFERENCE_SEQUENCE={reference.base} \
      OUTPUT=/io/prefix \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectQualityYieldMetrics \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE
      
    ls /io/
    cp /io/prefix.alignment_summary_metrics {j.out_alignment_summary_metrics}
    cp /io/prefix.base_distribution_by_cycle_metrics {j.out_base_distribution_by_cycle_metrics}
    cp /io/prefix.insert_size_metrics {j.out_insert_size_metrics}
    cp /io/prefix.quality_by_cycle_metrics {j.out_quality_by_cycle_metrics}
    cp /io/prefix.quality_yield_metrics {j.out_quality_yield_metrics}
    """

    j.command(wrap_command(cmd))
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
    out_hs_metrics_path: Path | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run picard CollectHsMetrics metrics.
    Based on https://github.com/broadinstitute/warp/blob/master/tasks/broad/Qc.wdl#L528
    """
    if can_reuse(out_hs_metrics_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'picard'}
    j = b.new_job('Picard CollectHsMetrics', job_attrs)
    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'exome'
    res = STANDARD.request_resources(mem_gb=7)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)
    targets_interval_file = b.read_input(
        str(reference_path('broad/exome_targets_interval_list'))
    )
    bait_interval_file = b.read_input(
        str(reference_path('broad/exome_bait_interval_list'))
    )

    cmd = f"""\
    grep -v 
    
    picard -Xms5g -Xmx{res.get_java_mem_mb()}m \
      CollectHsMetrics \
      INPUT={cram.cram} \
      REFERENCE_SEQUENCE={reference.base} \
      VALIDATION_STRINGENCY=SILENT \
      TARGET_INTERVALS={targets_interval_file} \
      BAIT_INTERVALS={bait_interval_file} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY \
      OUTPUT={j.out_hs_metrics}
    """

    j.command(wrap_command(cmd))
    b.write_output(j.out_hs_metrics, str(out_hs_metrics_path))
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

    job_attrs = (job_attrs or {}) | {'tool': 'picard'}
    j = b.new_job('Picard CollectWgsMetrics', job_attrs)

    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'genome'
    res = STANDARD.request_resources(ncpu=4)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)
    interval_file = b.read_input(
        str(reference_path('broad/genome_coverage_interval_list'))
    )

    cmd = f"""\
    picard -Xms2000m -Xmx{res.get_java_mem_mb()}m \
      CollectWgsMetrics \
      INPUT={cram.cram} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE={reference.base} \
      INTERVALS={interval_file} \
      OUTPUT={j.out_csv} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH={read_length}
    """

    j.command(wrap_command(cmd))
    b.write_output(j.out_csv, str(out_picard_wgs_metrics_path))
    return j
