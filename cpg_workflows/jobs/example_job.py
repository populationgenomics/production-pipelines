import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.hail_batch import command
from cpg_utils.hail_batch import image_path

def example_job(b: hb.Batch, output_path: Path, input_path: Path, job_attrs: dict | None = None) -> Job:
    j = b.new_job('copy', (job_attrs or {}) | {'tool': 'copying'})
    j.image(image_path('cpg_workflows'))
    input_file = b.read_input(str(input_path))
    cmd = f"""\
        echo {input_file} > {j.out}
    """
    j.command(command(cmd))
    b.write_output(j.out, str(output_path))
    return j
