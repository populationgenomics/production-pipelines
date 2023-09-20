"""
Convert BAM to CRAM.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import ResourceGroup
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command, image_path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD, HIGHMEM
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
    BamPath,
    CramPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
import re


def bam_to_cram(
    b: hb.Batch,
    input_bam: ResourceGroup,
    output_cram: CramPath | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceGroup] | tuple[None, CramPath]:
    """
    Convert a BAM file to a CRAM file.
    """
    # Don't run if output files exist and can be reused
    if output_cram and can_reuse(output_cram, overwrite):
        return None, output_cram
    
    assert isinstance(input_bam, ResourceGroup)

    job_name = 'bam_to_cram'
    if extra_label:
        job_name += f' {extra_label}'

    convert_tool = 'samtools_view'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=convert_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Get fasta file
    fasta_path = str(get_config()['references']['fasta'])
    fasta = b.read_input_group(
        fasta=fasta_path,
        fasta_fai=f'{fasta_path}.fai',
    )
    
    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    j.declare_resource_group(
        sorted_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )

    cmd = f'samtools view -@ {res.get_nthreads() - 1} -T {fasta.fasta} -C {input_bam.bam} | tee {j.sorted_cram["cram"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_cram["cram.crai"]}'
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_cram:
        output_bam_path = to_path(output_cram.path)
        b.write_output(j.sorted_cram, str(output_bam_path.with_suffix('')))

    return j, j.sorted_cram


def cram_to_bam(
    b: hb.Batch,
    input_cram: CramPath,
    output_bam: BamPath | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, ResourceGroup]:
    """
    Convert a CRAM file to a BAM file.
    """
    
    assert isinstance(input_cram, CramPath)

    cram_rg = b.read_input_group(
        cram=str(input_cram.path),
        crai=str(input_cram.index_path),
    )

    job_name = 'cram_to_bam'
    if extra_label:
        job_name += f' {extra_label}'

    convert_tool = 'samtools_view_cram_to_bam'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=convert_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Get fasta file
    # fasta_path = str(get_config()['references']['fasta'])
    # fasta = b.read_input_group(
    #     fasta=fasta_path,
    #     fasta_fai=f'{fasta_path}.fai',
    # )
    
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
        }
    )

    cmd = f'samtools view -@ {res.get_nthreads() - 1} -b {cram_rg.cram} | tee {j.sorted_bam["bam"]} | samtools index -@ {res.get_nthreads() - 1} - {j.sorted_bam["bam.bai"]}'
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_bam:
        output_bam_path = to_path(output_bam.path)
        b.write_output(j.sorted_bam, str(output_bam_path.with_suffix('')))

    return j, j.sorted_bam
