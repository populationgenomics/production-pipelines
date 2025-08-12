"""
Jobs that say hello.
"""

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import command, get_batch
from cpg_workflows.filetypes import CramPath


def hello_job(
    cram_path: CramPath,
    job_attrs: dict[str, str],
    output_base_path: Path,
) -> list[Job]:
    j = get_batch().new_job('Say hello', job_attrs)

    cmd = f"""
    echo Hello > {j.greet}
    """

    j.command(command(cmd))
    j.image(image_path('samtools'))
    get_batch().write_output(j.greet, str(output_base_path) + '.hello')

    j2 = get_batch().new_job('No-op', job_attrs)
    j2.command(command('echo hi there'))
    j2.image(image_path('samtools'))

    return [j, j2]
