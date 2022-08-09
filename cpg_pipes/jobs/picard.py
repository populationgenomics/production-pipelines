"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import HIGHMEM, STANDARD
from cpg_pipes.utils import can_reuse


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
) -> Job:
    """
    Make job that runs Picard CollectVariantCallingMetrics.
    """
    job_attrs = (job_attrs or {}) | dict(tool='picard_CollectVariantCallingMetrics')
    j = b.new_job('CollectVariantCallingMetrics', job_attrs)
    j.image(image_path('picard'))
    res = STANDARD.set_resources(j, storage_gb=20, mem_gb=3)
    reference = fasta_res_group(b)
    dbsnp_vcf = b.read_input_group(
        base=str(reference_path('broad/dbsnp_vcf')),
        index=str(reference_path('broad/dbsnp_vcf_index')),
    )
    intervals_file = b.read_input(
        str(reference_path('broad/genome_evaluation_interval_lists'))
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
