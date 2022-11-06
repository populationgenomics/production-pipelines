"""
Create Hail Batch jobs for VerifyBAMID.
"""
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from cpg_utils import Path
from cpg_utils.hail_batch import command
from cpg_workflows.utils import can_reuse
from cpg_workflows.filetypes import CramPath
from cpg_workflows.resources import STANDARD, storage_for_cram_qc_job

from hailtop.batch.job import Job


def verifybamid(
    b,
    cram_path: CramPath,
    out_verify_bamid_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run `VerifyBamID` contamination checks.
    Based on https://github.com/broadinstitute/warp/blob/57edec5591182d120b7d288b4b409e92a6539871/tasks/broad/BamProcessing.wdl#L395

    Creates a *.selfSM file, a TSV file with 2 rows, 19 columns.
    First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    """
    if can_reuse(out_verify_bamid_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'VerifyBamID'}
    j = b.new_job('VerifyBamID', job_attrs)
    j.image(image_path('verifybamid'))
    res = STANDARD.request_resources(ncpu=4)
    res.attach_disk_storage_gb = storage_for_cram_qc_job()
    res.set_to_job(j)
    reference = fasta_res_group(b)
    sequencing_type = get_config()['workflow']['sequencing_type']
    contam_ud = reference_path(f'broad/{sequencing_type}_contam_ud')
    contam_bed = reference_path(f'broad/{sequencing_type}_contam_bed')
    contam_mu = reference_path(f'broad/{sequencing_type}_contam_mu')
    if sequencing_type == 'exome':
        extra_opts = '--max-depth 1000'
    else:
        extra_opts = ''

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI
        
    /root/micromamba/share/verifybamid2-2.0.1-8/VerifyBamID \
    --NumThread {res.get_nthreads()} \
    --Verbose \
    --NumPC 4 \
    --Output OUTPUT \
    --BamFile $CRAM \
    --Reference {reference.base} \
    --UDPath {b.read_input(str(contam_ud))} \
    --MeanPath {b.read_input(str(contam_mu))} \
    --BedPath {b.read_input(str(contam_bed))} \
    {extra_opts}
    1>/dev/null

    cp OUTPUT.selfSM {j.out_selfsm}
    """

    j.command(
        command(
            cmd,
            define_retry_function=True,
        )
    )
    b.write_output(j.out_selfsm, str(out_verify_bamid_path))
    return j
