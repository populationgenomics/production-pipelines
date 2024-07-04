"""
Mark duplicates in BAM files using sambamba markdup
"""

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils.config import image_path
from cpg_utils.hail_batch import Batch, command
from cpg_workflows.filetypes import BamPath
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import Path


class Markdup:
    """
    Construct a sambamba markdup command for marking duplicates
    """

    def __init__(self, input_bam: BamPath | Path | str, output_bam: Path | str, nthreads: int = 8):
        self.command = [
            'sambamba markdup',
            '-t',
            str(nthreads),
            '--tmpdir=$BATCH_TMPDIR/markdup',
            str(input_bam),
            str(output_bam),
        ]

    def __str__(self) -> str:
        return 'mkdir -p $BATCH_TMPDIR/markdup && ' + ' '.join(self.command)

    def __repr__(self) -> str:
        return str(self)


def markdup(
    b: Batch,
    input_bam: ResourceGroup,
    job_attrs: dict | None = None,
    extra_label: str | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job | None, ResourceGroup]:
    """
    Takes an input BAM file and creates a job to mark duplicates with sambamba markdup.
    """

    assert isinstance(input_bam, ResourceGroup)

    base_job_name = 'sambamba_markdup'
    if extra_label:
        base_job_name += f' {extra_label}'

    tool = 'sambamba'

    j_name = base_job_name
    j_attrs = (job_attrs or {}) | dict(label=base_job_name, tool=tool)
    j = b.new_job(j_name, j_attrs)
    j.image(image_path('sambamba'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(j, ncpu=nthreads, storage_gb=50)

    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )
    assert isinstance(j.output_bam, ResourceGroup)

    cmd = Markdup(
        input_bam=input_bam.bam,
        output_bam=j.output_bam.bam,
        nthreads=res.get_nthreads(),
    )
    j.command(command(str(cmd), monitor_space=True))

    return j, j.output_bam
