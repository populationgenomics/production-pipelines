"""
Create Hail Batch jobs for samtools.
"""

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import command, fasta_res_group
from cpg_workflows.filetypes import CramPath
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse


def samtools_stats(
    b,
    cram_path: CramPath,
    out_samtools_stats_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run `samtools stats` for alignment QC.
    """
    if can_reuse(out_samtools_stats_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'samtools'}
    j = b.new_job('samtools stats', job_attrs)

    j.image(image_path('samtools'))
    res = STANDARD.set_resources(j, fraction=1)
    reference = fasta_res_group(b)

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    samtools stats \\
    -@{res.get_nthreads() - 1} \\
    --reference {reference.base} \\
    $CRAM > {j.output_stats}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_stats, str(out_samtools_stats_path))
    return j
