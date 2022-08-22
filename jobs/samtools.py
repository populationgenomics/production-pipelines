"""
Create Hail Batch jobs for samtools.
"""
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.hail_batch import image_path, fasta_res_group, command
from cpg_utils.flows.resources import STANDARD
from cpg_utils.flows.filetypes import CramPath
from cpg_utils.flows.utils import can_reuse


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
    cram = cram_path.resource_group(b)
    reference = fasta_res_group(b)

    cmd = f"""\
    samtools stats \\
    -@{res.get_nthreads() - 1} \\
    --reference {reference.base} \\
    {cram.cram} > {j.output_stats}
    """

    j.command(command(cmd))
    b.write_output(j.output_stats, str(out_samtools_stats_path))
    return j
