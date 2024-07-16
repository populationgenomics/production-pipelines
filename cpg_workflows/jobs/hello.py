"""
Test.
"""

from collections.abc import Iterable

from hailtop.batch.job import BashJob, Job
from hailtop.batch.resource import JobResourceFile, Resource, ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import command, fasta_res_group, get_batch, query_command
from cpg_workflows.filetypes import CramPath
from cpg_workflows.query_modules import seqr_loader, seqr_loader_cnv
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.scripts import upgrade_ped_with_inferred
from cpg_workflows.utils import can_reuse, chunks


def per_seqgroup(
    b,
    cram_path: CramPath,
    job_attrs: dict[str, str],
    outputs: dict[str, Path],
) -> list[Job]:
    j = b.new_bash_job('SeqGroup job', job_attrs)

    print(f'{type(cram_path.path)=}')
    cram = b.read_input(cram_path.path)

    cmd = f"""
    echo {cram}
    echo "Bump statistics (whatever they are) from {cram}" > {j.bumps}
    echo "Count statistics from {cram}" > {j.counts}
    """

    j.command(command(cmd))
    b.write_output(j.bumps, str(outputs['bump']))
    b.write_output(j.counts, str(outputs['count']))
    return [j]


def per_cohort(
    b,
    job_attrs: dict[str, str],
    output_path: Path,
) -> list[Job]:
    j = b.new_bash_job('Cohort job', job_attrs)

    cmd = f"""
    echo Hello {output_path} stuff
    """

    j.command(command(cmd))

    return [j]
