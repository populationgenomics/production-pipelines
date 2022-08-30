"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from cpg_utils.hail_batch import command
from cpg_utils.workflows.resources import HIGHMEM, STANDARD, storage_for_cram_qc_job
from cpg_utils.workflows.filetypes import CramPath
from cpg_utils.workflows.utils import can_reuse, exists


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
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)

    cmd = f"""\
    picard -Xmx{res.get_java_mem_mb()}m \
      CollectMultipleMetrics \
      INPUT={cram.cram} \
      REFERENCE_SEQUENCE={reference.base} \
      OUTPUT=$BATCH_TMPDIR/prefix \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectQualityYieldMetrics \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE
      
    ls $BATCH_TMPDIR/
    cp $BATCH_TMPDIR/prefix.alignment_summary_metrics {j.out_alignment_summary_metrics}
    cp $BATCH_TMPDIR/prefix.base_distribution_by_cycle_metrics {j.out_base_distribution_by_cycle_metrics}
    cp $BATCH_TMPDIR/prefix.insert_size_metrics {j.out_insert_size_metrics}
    cp $BATCH_TMPDIR/prefix.quality_by_cycle_metrics {j.out_quality_by_cycle_metrics}
    cp $BATCH_TMPDIR/prefix.quality_yield_metrics {j.out_quality_yield_metrics}
    """

    j.command(command(cmd))
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

    job_attrs = (job_attrs or {}) | {'tool': 'picard_CollectHsMetrics'}
    j = b.new_job('Picard CollectHsMetrics', job_attrs)
    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'exome'
    res = STANDARD.request_resources(ncpu=2)
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
    
    picard -Xmx{res.get_java_mem_mb()}m \
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

    j.command(command(cmd))
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

    job_attrs = (job_attrs or {}) | {'tool': 'picard_CollectWgsMetrics'}
    j = b.new_job('Picard CollectWgsMetrics', job_attrs)

    j.image(image_path('picard'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    assert sequencing_type == 'genome'
    res = STANDARD.request_resources(ncpu=2)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)
    interval_file = b.read_input(
        str(reference_path('broad/genome_coverage_interval_list'))
    )

    cmd = f"""\
    picard -Xmx{res.get_java_mem_mb()}m \
      CollectWgsMetrics \
      INPUT={cram.cram} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE={reference.base} \
      INTERVALS={interval_file} \
      OUTPUT={j.out_csv} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH={read_length}
    """

    j.command(command(cmd))
    b.write_output(j.out_csv, str(out_picard_wgs_metrics_path))
    return j
