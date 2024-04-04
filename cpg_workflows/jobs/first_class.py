"""
Create Hail Batch jobs to run STRipy
"""

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from cpg_workflows.resources import STANDARD
from cpg_workflows.targets import SequencingGroup


def first_class(
    b,
    sequencing_group: SequencingGroup,
    out_path: Path,
    job_attrs: dict | None = None,
) -> Job | None:
    """
    Run First Class Test
    """

    job_attrs = job_attrs or {}
    j = b.new_job('First Class File Check', job_attrs)

    cmd = f"""\
    echo "Running first class test for {sequencing_group.id}" > {j.ofile}
    """

    j.command(command(cmd))
    b.write_output(j.ofile, str(out_path))
    return j
