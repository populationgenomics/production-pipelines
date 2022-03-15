"""
Create Hail Batch jobs for alignment QC.
"""
import logging
from cloudpathlib import CloudPath
from hailtop.batch.job import Job

from cpg_pipes import images, buckets, ref_data
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.pipeline.analysis import CramPath

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def samtools_stats(
    b,
    cram_path: CramPath,
    sample_name: str,
    output_path: CloudPath | None = None,
    dataset_name: str | None = None,
    overwrite: bool = True,
) -> Job:
    """
    Run `samtools stats` for alignment QC.
    """
    jname = 'samtools stats'
    j = b.new_job(jname, dict(sample=sample_name, dataset=dataset_name))
    if not output_path:
        output_path = cram_path.path.with_suffix('.stats')
    if buckets.can_reuse(cram_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.SAMTOOLS_PICARD_IMAGE)
    job_resource = STANDARD.set_resources(j, storage_gb=60)
    cram = cram_path.resource_group(b)
    reference = b.read_input_group(**ref_data.REF_D)

    cmd = f"""\
    samtools stats \\
    -@{job_resource.get_nthreads() - 1} \\
    --reference {reference.base} \\
    {cram.cram} > {j.output_stats}
    """

    j.command(wrap_command(cmd))
    b.write_output(j.output_stats, str(output_path))
    return j


def verify_bamid(
    b,
    cram_path: CramPath,
    sample_name: str,
    output_path: CloudPath | None = None,
    dataset_name: str | None = None,
    overwrite: bool = True,
) -> Job:
    """
    Run `VerifyBamID` contamination checks. 
    Based on https://github.com/broadinstitute/warp/blob/57edec5591182d120b7d288b4b409e92a6539871/tasks/broad/BamProcessing.wdl#L395
    """
    jname = 'VerifyBamID'
    j = b.new_job(jname, dict(sample=sample_name, dataset=dataset_name))
    if not output_path:
        output_path = cram_path.path.with_suffix('.selfSM')
    if buckets.can_reuse(cram_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.VERIFY_BAMID)
    STANDARD.set_resources(j, storage_gb=60)
    cram = cram_path.resource_group(b)
    reference = b.read_input_group(**ref_data.REF_D)
    cont_ref_d = dict(
        ud=CloudPath(ref_data.CONTAM_BUCKET) / '1000g.phase3.100k.b38.vcf.gz.dat.UD',
        bed=CloudPath(ref_data.CONTAM_BUCKET) / '1000g.phase3.100k.b38.vcf.gz.dat.bed',
        mu=CloudPath(ref_data.CONTAM_BUCKET) / '1000g.phase3.100k.b38.vcf.gz.dat.mu',
    )
    res_group = b.read_input_group(**{k: str(v) for k, v in cont_ref_d.items()})
    
    output_prefix = sample_name
    cmd = f"""\
    # creates a *.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output {output_prefix} \
    --BamFile {cram.cram} \
    --Reference {reference.base} \
    --UDPath {res_group.ud} \
    --MeanPath {res_group.mu} \
    --BedPath {res_group.bed} \
    1>/dev/null
    
    cp {output_prefix}.selfSM {j.out_selfsm}
    """

    j.command(wrap_command(cmd))
    b.write_output(j.out_selfsm, str(output_path))
    return j


def picard_wgs_metrics(
    b,
    cram_path: CramPath,
    sample_name: str,
    output_path: CloudPath | None = None,
    dataset_name: str | None = None,
    overwrite: bool = True,
    read_length: int = 250,
) -> Job:
    """
    Run PicardTools CollectWgsMetrics metrics. 
    Based on https://github.com/broadinstitute/warp/blob/e1ac6718efd7475ca373b7988f81e54efab608b4/tasks/broad/Qc.wdl#L444
    """
    jname = 'Picard CollectWgsMetrics'
    j = b.new_job(jname, dict(sample=sample_name, dataset=dataset_name))
    if not output_path:
        output_path = cram_path.path.with_suffix('.csv')
    if buckets.can_reuse(cram_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.SAMTOOLS_PICARD_IMAGE)
    res = STANDARD.set_resources(j, storage_gb=60)
    cram = cram_path.resource_group(b)
    reference = b.read_input_group(**ref_data.REF_D)
    interval_file = b.read_input(ref_data.WGS_COVERAGE_INTERVAL_LIST)

    cmd = f"""\
    picard -Xms2000m -Xmx{res.get_java_mem_mb()}g \
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
    b.write_output(j.out_csv, str(output_path))
    return j
