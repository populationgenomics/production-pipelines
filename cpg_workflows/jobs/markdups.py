"""
Mark duplicates in BAM files using sambamba markdup
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import ResourceFile, ResourceGroup
from cpg_utils.hail_batch import command, image_path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse, Path, to_path
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
from dataclasses import dataclass
from enum import Enum


class Markdup:
    """
    Construct a sambamba markdup command for marking duplicates
    """

    def __init__(
        self,
        input_bam: BamPath | Path | str,
        output_bam: Path | str,
        nthreads: int = 8,
    ):
        self.command = [
            'sambamba markdup',
            '-t', str(nthreads),
            f'--tmpdir=$BATCH_TMPDIR/markdup',
            str(input_bam),
            str(output_bam),
        ]

    def __str__(self) -> str:
        return 'mkdir -p $BATCH_TMPDIR/markdup && ' + ' '.join(self.command)
    
    def __repr__(self) -> str:
        return str(self)


def markdup(
    b: hb.Batch,
    input_bam: ResourceGroup,
    output_bam: Path | str,
    job_attrs: dict | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> Job | None:
    """
    Takes an input BAM file and creates a job to mark duplicates with sambamba markdup.
    """
    # Don't run if output files exist and can be reused
    if output_bam and can_reuse(output_bam, overwrite):
        return None
    
    assert isinstance(input_bam, ResourceGroup)
    
    base_job_name = 'sambamba_markdup'
    if extra_label:
        base_job_name += f' {extra_label}'
    
    tool = 'sambamba'

    j_name = base_job_name
    j_attrs = (job_attrs or {}) | dict(label=base_job_name, tool=tool)
    j = b.new_job(j_name, j_attrs)
    # j.image(image_path('sambamba'))
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/sambamba:1.0.1')
    
    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    j.declare_resource_group(output_bam={
        'bam': '{root}.bam',
        'bam.bai': '{root}.bam.bai',
    })
    assert isinstance(j.output_bam, ResourceGroup)

    cmd = Markdup(
        input_bam=input_bam.bam,
        output_bam=j.output_bam.bam,
        nthreads=res.get_nthreads(),
    )
    j.command(command(str(cmd), monitor_space=True))

    # Write output to file
    if output_bam:
        output_bam_path = to_path(output_bam)
        b.write_output(j.output_bam, str(output_bam_path.with_suffix('')))

    return j
