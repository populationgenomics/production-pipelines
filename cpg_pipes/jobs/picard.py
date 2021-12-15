"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

from os.path import splitext, join, dirname
from typing import Optional

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import resources, utils
from cpg_pipes.hailbatch import wrap_command


def markdup(
    b: hb.Batch,
    sorted_bam: hb.ResourceFile,
    sample_name: str,
    project_name: Optional[str] = None,
    output_path: Optional[str] = None,
    overwrite: bool = True,
) -> Job:
    """
    Make job that runs Picard MarkDuplicates and converts the result to CRAM.
    """
    j = b.new_job('MarkDuplicates', dict(sample=sample_name, project=project_name))
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.SAMTOOLS_PICARD_IMAGE)
    ncpu = 4
    nthreads = ncpu * 2
    j.cpu(ncpu)  # ~3.75G per core, total 15G
    j.storage('200G')  # Allow 100G for input and 100G for output
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )
    fasta_reference = b.read_input_group(**resources.REF_D)

    cmd = f"""
    picard MarkDuplicates -Xms13G \\
    I={sorted_bam} O=/dev/stdout M={j.duplicate_metrics} \\
    TMP_DIR=$(dirname {j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate \\
    | samtools view -@{nthreads - 1} -T {fasta_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@{nthreads - 1} {j.output_cram.cram} {j.output_cram['cram.crai']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if output_path:
        b.write_output(j.output_cram, splitext(output_path)[0])
        b.write_output(
            j.duplicate_metrics,
            join(
                dirname(output_path),
                'duplicate-metrics',
                f'{sample_name}-duplicate-metrics.csv',
            ),
        )
    return j
