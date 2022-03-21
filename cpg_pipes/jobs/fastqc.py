"""
Jobs to run FastQC
"""
import hailtop.batch as hb
from hailtop.batch import ResourceFile
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.types import AlignmentInput, CramPath, FastqPair
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs import align
from cpg_pipes.refdata import RefData


def fastqc(
    b: hb.Batch, 
    output_html_path: Path,
    output_zip_path: Path,
    alignment_input: AlignmentInput,
    refs: RefData,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Adds FastQC jobs. If the input is a set of fqs, runs FastQC on each fq file.
    """
    def _fastqc_one(jname_, inp: ResourceFile):
        j = b.new_job(jname_, job_attrs)
        j.image('biocontainers/fastqc:v0.11.9_cv8')
        res = STANDARD.set_resources(j, ncpu=16)
        
        cmd = f"""\
        mkdir -p /io/batch/outdir
        fastqc -t{res.get_nthreads()} {inp} \\
        --outdir /io/batch/outdir
        ls /io/batch/outdir
        ln /io/batch/outdir/*_fastqc.html {j.out_html}
        ln /io/batch/outdir/*_fastqc.zip {j.out_zip}
        unzip /io/batch/outdir/*_fastqc.zip
        """
        j.command(wrap_command(cmd, monitor_space=True))
    
        b.write_output(j.out_html, str(output_html_path))
        b.write_output(j.out_zip, str(output_zip_path))
        return j
    
    jobs = []
    if isinstance(alignment_input, CramPath) and alignment_input.is_bam:
        bam = alignment_input.resource_group(b)
        j = _fastqc_one('FastQC', bam.bam)
        jobs.append(j)
        return jobs

    if isinstance(alignment_input, CramPath):
        extract_j = align.extract_fastq(
            b=b,
            cram=alignment_input.resource_group(b),
            refs=refs,
            job_attrs=job_attrs,
        )
        fastq_resouces = [FastqPair(extract_j.fq1, extract_j.fq2)]
    else:
        fastq_resouces = [pair.as_resources(b) for pair in alignment_input]

    for lane_i, pair in enumerate(fastq_resouces):
        jname = f'FastQC R1'
        if len(fastq_resouces) > 1:
            jname += f' lane={lane_i}'
        j = _fastqc_one(jname, pair.r1)
        jobs.append(j)
    return jobs
