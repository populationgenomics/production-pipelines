"""
Perform outlier gene expression analysis with Outrider.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command, image_path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
from textwrap import dedent


class Outrider:
    """
    Construct an outrider command for performing outlier gene expression analysis.
    """

    def __init__(
            self,
            output_path: str | Path | None = None,
        ) -> None:
        self.output_path = output_path
        self.command = f"""
        Rscript --vanilla <<EOF
        library(OUTRIDER)
        df <- read.table("{self.output_path}", header=TRUE, sep="\\\\t")
        EOF
        """
        self.command = dedent(self.command).strip()

    def __str__(self):
        return self.command
    
    def __repr__(self):
        return self.__str__()


def outrider(
    b: hb.Batch,
    output_path: str | Path | None = None,
    sample_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
) -> Job:
    """
    Run Outrider.
    """

    # Create job
    job_name = f'outrider_{sample_name}' if sample_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='outrider')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('outrider'))

    # Create counting command
    outrider = Outrider(
        output_path=j.output,
    )

    cmd = str(outrider)
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_path:
        # NOTE: j.output is just a placeholder
        b.write_output(j.output, str(output_path))
    
    return j
