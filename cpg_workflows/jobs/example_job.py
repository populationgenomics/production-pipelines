"""
A job to copy a fastq file and rename it
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.hail_batch import command
from cpg_utils.hail_batch import image_path

def example_job(
    b: hb.Batch,
    input_path: Path,
    output_path: Path,
    job_attrs: dict | None = None
) -> Job:
    j = b.new_job('Rename fastq file', (job_attrs or {}) | {'tool': 'copying'})
    j.image(image_path('cpg_workflows'))

    # input_file
    input_file: str | hb.ResourceFile = b.read_input(str(input_path))

    cmd = (
        f'echo {input_file} > {j.out}'
    )
    j.command(command(cmd))

    b.write_output(j.out, str(output_path))
