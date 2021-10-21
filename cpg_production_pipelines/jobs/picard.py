from typing import Optional

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines.jobs import wrap_command
from cpg_production_pipelines import resources


def markdup(
    b: hb.Batch,
    sorted_bam: hb.ResourceFile,
    sample_name: str,
    project_name: str,
    out_cram_path: Optional[str] = None,
    overwrite: bool = True,
) -> Job:
    """
    Make job that runs Picard MarkDuplicates and converts the result to CRAM.
    """
    job_name = f'{sample_name}: MarkDuplicates'
    if project_name:
        job_name = f'{project_name}/{job_name}'
    j = b.new_job(job_name)
    j.image(resources.SAMTOOLS_PICARD_IMAGE)
    j.cpu(2)
    j.memory('highmem')
    j.storage('400G')  # Allow 200G for input and 200G for output
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )
    fasta_reference = b.read_input_group(**resources.REF_D)

    cmd = f"""
    picard MarkDuplicates -Xms7G \\
    I={sorted_bam} O=/dev/stdout M={j.duplicate_metrics} \\
    TMP_DIR=$(dirname {j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate \\
    | samtools view -@2 -T {fasta_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@2 {j.output_cram.cram} {j.output_cram['cram.crai']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    return j
