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
from os.path import basename


class Outrider:
    """
    Construct an outrider command for performing outlier gene expression analysis.
    """

    def __init__(
            self,
            input_counts: list[str | Path],
            output_path: str | Path | None = None,
        ) -> None:
        self.input_counts = input_counts
        assert isinstance(self.input_counts, list), f'input_counts must be a list, instead got {self.input_counts}'
        self.input_counts_r_str = ', '.join([f'"{str(f)}"' for f in self.input_counts])
        self.output_path = output_path
        self.command = f"""
        Rscript --vanilla <<EOF
        library(OUTRIDER)
        input_counts <- c({self.input_counts_r_str})
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
    input_counts: list[str | Path],
    output_path: str | Path | None = None,
    cohort_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
) -> Job:
    """
    Run Outrider.
    """

    # Localise input files
    assert all([isinstance(f, (str, Path)) for f in input_counts])
    infiles = {
        basename(str(f)).replace('.count', ''): b.read_input(str(f))
        for f in input_counts
    }

    # Create job
    job_name = f'outrider_{cohort_name}' if cohort_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='outrider')
    j = b.new_job(job_name, _job_attrs)
    # j.image(image_path('outrider'))
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/outrider:1.18.1')

    # Create counting command
    outrider = Outrider(
        input_counts=input_counts,
        output_path=j.output,
    )

    cmd = str(outrider)
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_path:
        # NOTE: j.output is just a placeholder
        b.write_output(j.output, str(output_path))
    
    return j
