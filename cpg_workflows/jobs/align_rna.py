"""
Align RNA-seq reads to the genome using STAR.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
    BamPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
import re


class STAR:
    """
    Construct a STAR command for aligning FASTQs.
    """

    def __init__(
            self,
            input_fastq_pair: FastqPair,
            sample_name: str,
            genome_dir: str | Path,
            nthreads: int,
            read_group: dict[str, str] | None = None,
            output_path: str | BamPath | Path | None = None,
            bamout: bool = True,
            sort: bool = True,
            stdout: bool = False,
    ):
        self.command = ['STAR']
        # Create read group line
        self.read_group = read_group or {}
        if 'ID' not in self.read_group:
            self.read_group['ID'] = sample_name
        if 'SM' not in self.read_group:
            self.read_group['SM'] = sample_name
        self.read_group_line = self._create_read_group_line()
        # Create outSAMtype and outStd strings
        if bamout and sort:
            self.outSAMtype = 'BAM SortedByCoordinate'
            self.outStd = 'BAM_SortedByCoordinate'
            self.output = 'Aligned.sortedByCoord.out.bam'
        elif bamout:
            self.outSAMtype = 'BAM Unsorted'
            self.outStd = 'BAM_Unsorted'
            self.output = 'Aligned.out.bam'
        else:
            self.outSAMtype = 'SAM'
            self.outStd = 'SAM'
            self.output = 'Aligned.out.sam'
        if not stdout:
            self.outStd = 'Log'
        # Get output extension
        self.output_extension = 'bam' if bamout else 'sam'
        # Create command
        self.command.extend([
            '--runThreadN', str(nthreads),
            '--genomeDir', str(genome_dir),
            '--outSAMtype', self.outSAMtype,
            '--outStd', self.outStd,
            '--outSAMattrRGline', self.read_group_line,
            '--readFilesIn', str(input_fastq_pair.r1), str(input_fastq_pair.r2),
        ])
        if output_path and not stdout:
            self.move_output_command = f' && mv {self.output} {output_path}'
        elif output_path and stdout:
            self.move_output_command = f' > {output_path}'
        else:
            self.move_output_command = ''

    def __str__(self):
        return ' '.join(self.command) + self.move_output_command
    
    def __repr__(self):
        return self.__str__()

    def _create_read_group_line(self) -> str:
        """
        Create a read group line for the STAR command.
        """
        read_group_line = ''
        for k, v in self.read_group.items():
            v_quoted = re.sub(r'(^.*\s.*$)', r'"\1"', v)
            read_group_line += f'{k}:{v_quoted} '
        return read_group_line.strip()
    

def align(
    b: hb.Batch,
    fastq_pairs: FastqPairs,
    output_bam: BamPath,
    sample_name: str,
    genome_dir: str | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
):
    """
    Align (potentially multiple) FASTQ pairs using STAR,
    merge the resulting BAMs (if necessary),
    then sort and index the resulting BAM.
    """
    # Don't run if the output exists and can be reused
    if output_bam and can_reuse(output_bam, overwrite):
        return None
    
    if not isinstance(fastq_pairs, FastqPairs):
        raise TypeError(f'fastq_pairs must be a FastqPairs object, not {type(fastq_pairs)}')
    if len(fastq_pairs) == 0:
        raise ValueError('fastq_pairs must contain at least one FastqPair')
    
    merge = len(fastq_pairs) > 1

    jobs = []
    aligned_bams = []
    job_idx = 1
    for fq_pair in fastq_pairs:
        if not isinstance(fq_pair, FastqPair):
            raise TypeError(f'fastq_pairs must contain FastqPair objects, not {type(fq_pair)}')
        label = f'{extra_label} {job_idx}' if extra_label else f'{job_idx}'
        j = align_fq_pair(
            b=b,
            fastq_pair=fq_pair,
            sample_name=sample_name,
            genome_dir=genome_dir,
            extra_label=label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        aligned_bams.append(j.output_bam)
        job_idx += 1
    
    if merge:
        j = merge_bams(
            b=b,
            input_bams=aligned_bams,
            extra_label=extra_label,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(j)
        aligned_bam = j.merged_bam
    else:
        aligned_bam = aligned_bams[0]

    j = sort_index_bam(
        b=b,
        input_bam=aligned_bam,
        extra_label=extra_label,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )
    jobs.append(j)

    if output_bam:
        b.write_output(j.sorted_bam, str(output_bam.path))

    return jobs


def align_fq_pair(
    b: hb.Batch,
    fastq_pair: FastqPair,
    sample_name: str,
    genome_dir: str | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> Job:
    """
    Takes an input FastqPair object, and creates a job to align it using STAR.
    """
    job_name = 'align_rna'
    if extra_label:
        job_name += f' {extra_label}'

    # If number of threads is not specified, use a whole instance
    nthreads = requested_nthreads or STANDARD.max_threads()

    align_tool = 'STAR'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=align_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    star = STAR(
        input_fastq_pair=fastq_pair,
        sample_name=sample_name,
        genome_dir=genome_dir,
        nthreads=(nthreads - 1),
        output_path=j.output_bam,
        bamout=True,
        sort=True,
        stdout=False,
    )
    cmd = str(star)
    j.command(command(cmd, monitor_space=True))
    return j
    

def merge_bams(
    b: hb.Batch,
    input_bams: list[str | BamPath | Path],
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
) -> Job:
    """
    Merge a list of BAM files into a single BAM file.
    """
    job_name = 'merge_bams'
    if extra_label:
        job_name += f' {extra_label}'
    
    # If number of threads is not specified, use a whole instance
    nthreads = requested_nthreads or STANDARD.max_threads()

    merge_tool = 'samtools'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=merge_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    cmd = f'samtools merge -@ {nthreads - 1} -o {j.merged_bam} {" ".join([str(b) for b in input_bams])}'
    j.command(command(cmd, monitor_space=True))
    return j


def sort_index_bam(
    b: hb.Batch,
    input_bam: str | BamPath | Path,
    extra_label: str | None = None,
    job_attrs: dict | None = None,
    requested_nthreads: int | None = None,
):
    """
    Sort and index a BAM file.
    """
    job_name = 'sort_index_bam'
    if extra_label:
        job_name += f' {extra_label}'
    
    # If number of threads is not specified, use a whole instance
    nthreads = requested_nthreads or STANDARD.max_threads()

    sort_tool = 'samtools'
    j_attrs = (job_attrs or {}) | dict(label=job_name, tool=sort_tool)
    j = b.new_job(name=job_name, attributes=j_attrs)
    cmd = f'samtools sort -@ {nthreads - 1} {input_bam} | tee {j.sorted_bam} | samtools index -@ {nthreads - 1} - {j.sorted_bam_idx}'
    j.command(command(cmd, monitor_space=True))
    return j
