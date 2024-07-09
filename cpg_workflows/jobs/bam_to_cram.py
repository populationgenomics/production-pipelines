"""
Convert BAM to CRAM.
"""

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import Batch, command
from cpg_workflows.resources import STANDARD


def bam_to_cram(
    b: Batch,
    input_bam: ResourceGroup,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
    reference_fasta_path: Path | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a BAM file to a CRAM file.
    """

    assert isinstance(input_bam, ResourceGroup)

    job_name = 'bam_to_cram'
    if extra_label:
        job_name += f' {extra_label}'

    convert_tool = 'samtools_view'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=convert_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Get fasta file
    fasta = b.read_input_group(
        fasta=reference_fasta_path,
        fasta_fai=f'{reference_fasta_path}.fai',
    )

    # Set resource requirements
    nthreads = requested_nthreads or 8
    # TODO: make storage configurable
    res = STANDARD.set_resources(j, ncpu=nthreads, storage_gb=50)

    j.declare_resource_group(
        sorted_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )

    cmd = f'samtools view -@ {res.get_nthreads() - 1} -T {fasta.fasta} -C {input_bam.bam} | tee {j.sorted_cram["cram"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_cram["cram.crai"]}'
    j.command(command(cmd, monitor_space=True))

    return j, j.sorted_cram


def cram_to_bam(
    b: Batch,
    input_cram: ResourceGroup,
    output_bam: Path | None = None,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a CRAM file to a BAM file.
    """

    assert isinstance(input_cram, ResourceGroup)

    job_name = 'cram_to_bam'
    if extra_label:
        job_name += f' {extra_label}'

    convert_tool = 'samtools_view_cram_to_bam'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=convert_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    j.declare_resource_group(
        sorted_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        },
    )

    cmd = f'samtools view -@ {res.get_nthreads() - 1} -b {input_cram.cram} | tee {j.sorted_bam["bam"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_bam["bam.bai"]}'
    j.command(command(cmd, monitor_space=True))

    # Write BAM if requested
    if output_bam:
        b.write_output(j.sorted_bam, str(output_bam.with_suffix('')))

    return j, j.sorted_bam
