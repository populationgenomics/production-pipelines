"""
Create Hail Batch jobs for alignment QC.
"""
import logging

from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.types import CramPath
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD

logger = logging.getLogger(__file__)


def samtools_stats(
    b,
    cram_path: CramPath,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Run `samtools stats` for alignment QC.
    """
    jname = 'samtools stats'
    j = b.new_job(jname, job_attrs)
    if not output_path:
        output_path = cram_path.path.with_suffix('.stats')
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(image_path('samtools'))
    job_resource = STANDARD.set_resources(j, storage_gb=60)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)

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
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Run `VerifyBamID` contamination checks.
    Based on https://github.com/broadinstitute/warp/blob/57edec5591182d120b7d288b4b409e92a6539871/tasks/broad/BamProcessing.wdl#L395
    """
    jname = 'VerifyBamID'
    j = b.new_job(jname, job_attrs)
    if not output_path:
        output_path = cram_path.path.with_suffix('.selfSM')
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(image_path('verify-bam-id'))
    STANDARD.set_resources(j, storage_gb=60)
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)
    contam_ud = reference_path('broad/contam_ud')
    contam_bed = reference_path('broad/contam_bed')
    contam_mu = reference_path('broad/contam_mu')

    cmd = f"""\
    # creates a *.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output OUTPUT \
    --BamFile {cram.cram} \
    --Reference {reference.base} \
    --UDPath {b.read_input(contam_ud)} \
    --MeanPath {b.read_input(contam_mu)} \
    --BedPath {b.read_input(contam_bed)} \
    1>/dev/null

    cp OUTPUT.selfSM {j.out_selfsm}
    """

    j.command(wrap_command(cmd))
    b.write_output(j.out_selfsm, str(output_path))
    return j


def picard_wgs_metrics(
    b,
    cram_path: CramPath,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
    read_length: int = 250,
) -> Job:
    """
    Run PicardTools CollectWgsMetrics metrics.
    Based on https://github.com/broadinstitute/warp/blob/e1ac6718efd7475ca373b7988f81e54efab608b4/tasks/broad/Qc.wdl#L444
    """
    jname = 'Picard CollectWgsMetrics'
    j = b.new_job(jname, job_attrs)
    if not output_path:
        output_path = cram_path.path.with_suffix('.csv')
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(image_path('picard'))
    res = STANDARD.set_resources(j, storage_gb=60)
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
    b.write_output(j.out_csv, str(output_path))
    return j
