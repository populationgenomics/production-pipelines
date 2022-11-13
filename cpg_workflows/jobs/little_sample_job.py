"""
Sample Jobs.
"""
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.hail_batch import command
from cpg_utils.hail_batch import image_path


def little_sample_job(
    b: hb.Batch,
    output_path: Path,
    input_path: Path,
    subsample: bool = True,
    job_attrs: dict | None = None,
) -> Job:
    """
    Creates a copy of a file
    """

    j = b.new_job('TESTING WORKSHOP 22', (job_attrs or {}) | {'tool': 'silly_sample'})
    j.image(image_path('cpg_workflows'))

    cmd = ''
    input_file: str | hb.ResourceFile = b.read_input(str(input_path))

    cmd += f"""\
    echo {input_file} > {j.out}
    """
    j.command(command(cmd))
    b.write_output(j.out, str(output_path))
    return j
