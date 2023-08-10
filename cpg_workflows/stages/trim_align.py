"""
Align RNA-seq reads to the genome using STAR.
"""

import logging
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import ExpectedResultT
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroupStage,
)
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
    BamPath,
)
from cpg_workflows.jobs import trim
from cpg_workflows.jobs import align_rna
import re
from os.path import basename
from dataclasses import dataclass


def get_trim_inputs(sequencing_group: SequencingGroup) -> FastqPairs | None:
    """
    Get the input FASTQ file pairs for trimming
    """
    sequencing_type = get_config()['workflow']['sequencing_type']
    alignment_input = sequencing_group.alignment_input_by_seq_type.get(sequencing_type)
    if (
        not alignment_input or
        (get_config()['workflow'].get('check_inputs', True) and not alignment_input.exists()) or
        not isinstance(alignment_input, (FastqPair, FastqPairs))
    ):
        return None
    if isinstance(alignment_input, FastqPair):
        alignment_input = FastqPairs([alignment_input])
    return alignment_input


@dataclass
class InOutFastqPair:
    """
    Represents a single set of input and output paired FASTQ files
    """

    id: str
    input_pair: FastqPair
    output_pair: FastqPair


def get_input_output_pairs(sequencing_group: SequencingGroup) -> list[InOutFastqPair]:
    """
    Get the input FASTQ pairs, determine the trimmed output FASTQ file pair paths and
    output a list of InOutFastqPair objects
    """
    inputs = get_trim_inputs(sequencing_group)
    if not inputs or not isinstance(inputs, FastqPairs):
        return []
    prefix = sequencing_group.dataset.tmp_prefix() / 'trim'
    trim_suffix = '.trimmed.fastq.gz'
    input_output_pairs = []
    for i, pair in enumerate(inputs, 1):
        assert isinstance(pair, FastqPair)
        input_r1_bn = re.sub('.f(ast)?q.gz', '', basename(pair.r1))
        input_r2_bn = re.sub('.f(ast)?q.gz', '', basename(pair.r2))
        output_r1 = prefix / f'{input_r1_bn}{trim_suffix}'
        output_r2 = prefix / f'{input_r2_bn}{trim_suffix}'
        input_output_pairs.append(InOutFastqPair(
            str(i),  # ID
            pair,  # input FASTQ pair
            FastqPair(output_r1, output_r2),  # output FASTQ pair
        ))
        i += 1
    return input_output_pairs


@stage
class TrimAlignRNA(SequencingGroupStage):
    """
    Trim and align RNA-seq FASTQ reads with fastp and STAR
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Expect a pair of BAM and BAI files, one per set of input FASTQ files
        """
        return {
            suffix: sequencing_group.dataset.prefix() / 'bam' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('bam', 'bam'),
                ('bai', 'bam.bai'),
            ]
        }
    
    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to align the input FASTQ files to the genome using STAR
        """
        jobs = []

        # Run trim
        input_output_pairs = get_input_output_pairs(sequencing_group)
        for io_pair in input_output_pairs:
            try:
                j = trim.trim(
                    b=get_batch(),
                    sequencing_group=sequencing_group,
                    input_fq_pair=io_pair.input_pair,
                    output_fq_pair=io_pair.output_pair,
                    job_attrs=self.get_job_attrs(sequencing_group),
                    overwrite=sequencing_group.forced,
                    extra_label=f'fastq_pair_{io_pair.id}',
                )
                if j:
                    jobs.append(j)
            except trim.MissingFastqInputException:
                if get_config()['workflow'].get('skip_sgs_with_missing_input'):
                    logging.error(f'No FASTQ inputs, skipping sample {sequencing_group}')
                    sequencing_group.active = False
                    return self.make_outputs(sequencing_group, skipped=True)  # return empty output
                else:
                    return self.make_outputs(
                        target=sequencing_group, error_msg=f'No FASTQ input found'
                    )

        # Run alignment
        trimmed_fastq_pairs = FastqPairs([
            io_pair.output_pair
            for io_pair in input_output_pairs
        ])
        _exp_out = self.expected_outputs(sequencing_group)
        output_bam = BamPath(
            path=_exp_out['bam'],
            index_path=_exp_out['bai'],
        )
        try:
            jobs = align_rna.align(
                b=get_batch(),
                fastq_pairs=trimmed_fastq_pairs,
                output_bam=output_bam,
                sample_name=sequencing_group.id,
                genome_prefix=get_config()['references'].get('star_ref_dir'),
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
        except Exception as e:
            logging.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
        return self.make_outputs(
            sequencing_group,
            data=_exp_out,
            jobs=jobs
        )
