"""
Create Hail Batch jobs to run Picard tools (marking duplicates, QC).
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.providers.images import Images
from cpg_pipes.providers.refdata import RefData


def markdup(
    b: hb.Batch,
    sorted_bam: hb.ResourceFile,
    sample_name: str,
    refs: RefData,
    images: Images,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    qc_bucket: Path | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Make job that runs Picard MarkDuplicates and converts the result to CRAM.
    """
    job_attrs = (job_attrs or {}) | dict(tool='picard_MarkDuplicates')
    j = b.new_job('MarkDuplicates', job_attrs)
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.get('picard_samtools'))
    resource = STANDARD.set_resources(
        j, storage_gb=175
    )  # enough for input BAM and output CRAM
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )
    fasta_reference = refs.fasta_res_group(b)

    cmd = f"""
    picard MarkDuplicates -Xms13G \\
    I={sorted_bam} O=/dev/stdout M={j.duplicate_metrics} \\
    TMP_DIR=$(dirname {j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate \\
    | samtools view -@{resource.get_nthreads() - 1} -T {fasta_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@{resource.get_nthreads() - 1} {j.output_cram.cram} {j.output_cram['cram.crai']}
    """
    j.command(wrap_command(cmd, monitor_space=True))

    if output_path:
        if not qc_bucket:
            qc_bucket = output_path.parent

        b.write_output(j.output_cram, str(output_path.with_suffix('')))
        b.write_output(
            j.duplicate_metrics,
            str(
                qc_bucket / 'duplicate-metrics' / f'{sample_name}-duplicate-metrics.csv'
            ),
        )

    return j
