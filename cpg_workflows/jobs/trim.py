"""
Trim raw FASTQ reads using cutadapt
"""

from dataclasses import dataclass
from enum import Enum

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import Batch, command
from cpg_workflows.filetypes import FastqPair
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import SequencingGroup


class MissingFastqInputException(Exception):
    """Raise if alignment input is missing"""

    pass


class InvalidSequencingTypeException(Exception):
    """Raise if alignment type is not 'rna'"""

    pass


@dataclass
class AdapterSequence:
    """
    A class to represent an adapter sequence.
    """

    sequence: str
    name: str | None = None


@dataclass
class AdapterPair:
    """
    A class to represent a pair of adapter sequences.
    """

    r1: AdapterSequence
    r2: AdapterSequence


class AdapterPairs(Enum):
    """
    A class to represent a set of adapter pairs.
    """

    ILLUMINA_TRUSEQ = AdapterPair(
        r1=AdapterSequence(
            sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
            name='TruSeq Adapter Index 1',
        ),
        r2=AdapterSequence(
            sequence='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
            name='TruSeq Adapter Index 2',
        ),
    )


class Cutadapt:
    """
    Construct a cutadapt command for trimming FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        output_fastq_pair: FastqPair,
        adapter_type: str,
        paired: bool = True,
        min_length: int = 50,
        two_colour: bool = True,
        polyA: bool = False,
        quality_trim: int | None = None,  # TODO: determine default
    ):
        try:
            adapters: AdapterPair = AdapterPairs[adapter_type].value
        except AttributeError:
            raise ValueError(f'Invalid adapter type: {adapter_type}')
        self.command = [
            'cutadapt',
            *('-o', str(output_fastq_pair.r1)),
            *('-a', adapters.r1.sequence),
        ]
        if paired:
            self.command.extend(
                [
                    *('-p', str(output_fastq_pair.r2)),
                    *('-A', adapters.r2.sequence),
                ],
            )
        if quality_trim:
            if two_colour:
                self.command.append(f'--nextseq-trim={quality_trim}')
            else:
                self.command.append(f'-q {quality_trim}')
        if min_length:
            self.command.append(f'--minimum-length={min_length}')
        if polyA:
            self.command.append('--poly-a')
        self.command.append(str(input_fastq_pair.r1))
        if paired:
            self.command.append(str(input_fastq_pair.r2))

    def __str__(self) -> str:
        return ' '.join(self.command)

    def __repr__(self) -> str:
        return str(self)


class Fastp:
    """
    Construct a fastp command for trimming FASTQs.
    """

    def __init__(
        self,
        input_fastq_pair: FastqPair,
        output_fastq_pair: FastqPair,
        adapter_type: str,
        paired: bool = True,
        min_length: int = 50,
        nthreads: int = 3,
        polyG: bool = True,
        polyX: bool = False,
    ):
        try:
            adapters: AdapterPair = AdapterPairs[adapter_type].value
        except AttributeError:
            raise ValueError(f'Invalid adapter type: {adapter_type}')
        self.command = [
            'fastp',
            *('--in1', str(input_fastq_pair.r1)),
            *('--out1', str(output_fastq_pair.r1)),
            *('--length_required', str(min_length)),
            *('--adapter_sequence', adapters.r1.sequence),
            *('--thread', str(nthreads)),
        ]
        if paired:
            self.command.extend(
                [
                    *('--in2', str(input_fastq_pair.r2)),
                    *('--out2', str(output_fastq_pair.r2)),
                    *('--adapter_sequence_r2', adapters.r2.sequence),
                ],
            )
        if not polyG:
            self.command.append('--disable_trim_poly_g')
        if polyX:
            self.command.append('--trim_poly_x')

    def __str__(self) -> str:
        return ' '.join(self.command)

    def __repr__(self) -> str:
        return str(self)


def trim(
    b: Batch,
    sequencing_group: SequencingGroup,
    input_fq_pair: FastqPair,
    output_fq_pair: FastqPair | None = None,
    job_attrs: dict | None = None,
    extra_label: str | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> tuple[Job | None, FastqPair]:
    """
    Takes an input FastqPair object, and creates a job to trim the FASTQs using cutadapt.
    """
    # Don't run if all output files exist and can be reused
    if output_fq_pair and can_reuse(output_fq_pair.r1, overwrite) and can_reuse(output_fq_pair.r2, overwrite):
        return None, output_fq_pair.as_resources(b)

    base_job_name = 'TrimFastqs'
    if extra_label:
        base_job_name += f' {extra_label}'

    if not get_config()['workflow']['sequencing_type'] == 'transcriptome':
        raise InvalidSequencingTypeException(
            f"Invalid sequencing type '{get_config()['workflow']['sequencing_type']}'"
            + f" for job type '{base_job_name}'; sequencing type must be 'transcriptome'",
        )

    try:
        adapter_type = get_config()['trim']['adapter_type']
    except KeyError:
        raise ValueError('No adapter type specified in config file')

    trim_tool = 'fastp'

    trim_j_name = base_job_name
    trim_j_attrs = (job_attrs or {}) | dict(label=base_job_name, tool=trim_tool)
    trim_j = b.new_job(trim_j_name, trim_j_attrs)
    trim_j.image(image_path('fastp'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        trim_j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    fastq_pair = input_fq_pair.as_resources(b)

    trim_j.declare_resource_group(output_r1={'fastq.gz': '{root}.fastq.gz'})
    trim_j.declare_resource_group(output_r2={'fastq.gz': '{root}.fastq.gz'})
    assert isinstance(trim_j.output_r1, ResourceGroup)
    assert isinstance(trim_j.output_r2, ResourceGroup)
    out_fqs = FastqPair(
        r1=trim_j.output_r1['fastq.gz'],
        r2=trim_j.output_r2['fastq.gz'],
    )

    trim_cmd = Fastp(
        input_fastq_pair=fastq_pair,
        output_fastq_pair=out_fqs,
        adapter_type=adapter_type,
        paired=True,
        min_length=50,
        nthreads=res.get_nthreads(),
        polyG=True,
        polyX=True,
    )
    trim_j.command(command(str(trim_cmd), monitor_space=True))

    # Write output to file
    if output_fq_pair:
        b.write_output(out_fqs.r1, str(output_fq_pair.r1))
        b.write_output(out_fqs.r2, str(output_fq_pair.r2))

    return trim_j, out_fqs
