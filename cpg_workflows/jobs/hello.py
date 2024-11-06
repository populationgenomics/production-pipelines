"""
Jobs that say hello.
"""

from hailtop.batch.job import BashJob, Job
from hailtop.batch.resource import JobResourceFile, Resource, ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import command, fasta_res_group, get_batch, query_command
from cpg_workflows.filetypes import CramPath
from cpg_workflows.query_modules import seqr_loader
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.scripts import upgrade_ped_with_inferred
from cpg_workflows.utils import can_reuse, chunks


def hello_job(
    cram_path: CramPath,
    job_attrs: dict[str, str],
    output_base_path: Path,
) -> list[Job]:
    j = get_batch().new_job('Say hello', job_attrs | {'tool': 'echo'})

    cmd = f"""
    echo Hello > {j.greet}
    """

    j.command(command(cmd))
    get_batch().write_output(j.greet, str(output_base_path) + '.hello')
    return [j]
