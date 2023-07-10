"""
Trim raw FASTQ reads using cutadapt
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils.hail_batch import command
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)

class MissingFastqInputException(Exception):
    """Raise if alignment input is missing"""
    pass


class InvalidSequencingTypeException(Exception):
    """Raise if alignment type is not 'rna'"""
    pass


class Cutadapt:
    """
    Construct a cutadapt command for trimming FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        output_fastq_pair: FastqPair,
        adapterR13p: str | None = None,
        adapterR23p: str | None = None,
        adapterR15p: str | None = None,
        adapterR25p: str | None = None,
        quality_trim: int | None = None,  # TODO: determine default
        two_colour: bool = True,
        min_length: int | None = None,
        polyA: bool = False,
        poly: str | None = None,
        # TODO: add more parameters
    ):
        self.command = ['cutadapt', '-o', str(output_fastq_pair.r1), '-p', str(output_fastq_pair.r2)]
        adapter_args = []
        if adapterR13p:
            adapter_args.extend(['-a', adapterR13p])
        if adapterR23p:
            adapter_args.extend(['-A', adapterR23p])
        if adapterR15p:
            adapter_args.extend(['-g', adapterR15p])
        if adapterR25p:
            adapter_args.extend(['-G', adapterR25p])
        self.command.extend(adapter_args)
        if quality_trim:
            if two_colour:
                self.command.append(f'--nextseq-trim={quality_trim}')
            else:
                self.command.append(f'-q {quality_trim}')
        if min_length:
            self.command.append(f'--minimum-length={min_length}')
        if polyA:
            self.command.append('--poly-a')
        self.command.extend([str(input_fastq_pair.r1), str(input_fastq_pair.r2)])

    def __str__(self) -> str:
        return ' '.join(self.command)
    
    def __repr__(self) -> str:
        return str(self)


def trim(
    b: hb.Batch,
    sequencing_group: SequencingGroup,
    job_attrs: dict | None = None,
    input_fq_pair: FastqPair | None = None,
    output_fq_pair: FastqPair | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> Job:
    """
    Takes an input FastqPair object, and creates a job to trim the FASTQs using cutadapt.
    """
    # Don't run if all output files exist and can be reused
    if (
        output_fq_pair and
        can_reuse(output_fq_pair.r1, overwrite) and
        can_reuse(output_fq_pair.r2, overwrite)
    ):
        return None
    
    base_job_name = 'TrimFastqs'
    if extra_label:
        base_job_name += f' {extra_label}'
    
    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_threads()
    
    if not get_config()['workflow']['sequencing_type'] == 'rna':
        raise InvalidSequencingTypeException(
            f"Invalid sequencing type '{get_config()['workflow']['sequencing_type']}'" +
            f" for job type '{base_job_name}'; sequencing type must be 'rna'"
        )
    
    trim_tool = 'cutadapt'

    trim_j_name = base_job_name
    trim_j_attrs = (job_attrs or {}) | dict(label=base_job_name, tool=trim_tool)
    trim_j = b.new_job(trim_j_name, trim_j_attrs)
    trim_cmd = Cutadapt(  # TODO: add more arguments
        input_fastq_pair=input_fq_pair,
        output_fastq_pair=FastqPair(
            r1=trim_j.output_r1,
            r2=trim_j.output_r2,
        ),
    )
    trim_j.command(command(str(trim_cmd), monitor_space=True))

    # Write output to file
    if output_fq_pair:
        b.write_output(trim_j.output_r1, str(output_fq_pair.r1))
        b.write_output(trim_j.output_r2, str(output_fq_pair.r2))

    return trim_j