"""
Align RNA-seq reads to the genome using STAR.
"""

import logging
import re
from dataclasses import dataclass
from os.path import basename

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs
from cpg_workflows.jobs import align_rna, trim
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


def get_trim_inputs(sequencing_group: SequencingGroup) -> FastqPairs | None:
    """
    Get the input FASTQ file pairs for trimming
    """
    sequencing_type = get_config()['workflow']['sequencing_type']
    alignment_input = sequencing_group.alignment_input_by_seq_type.get(sequencing_type)
    if (
        not alignment_input
        or (get_config()['workflow'].get('check_inputs', True) and not alignment_input.exists())
        or not isinstance(alignment_input, (FastqPair, FastqPairs))
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
        input_r1_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r1)))
        input_r2_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r2)))
        output_r1 = prefix / f'{input_r1_bn}{trim_suffix}'
        output_r2 = prefix / f'{input_r2_bn}{trim_suffix}'
        input_output_pairs.append(
            InOutFastqPair(
                str(i),  # ID
                pair,  # input FASTQ pair
                FastqPair(output_r1, output_r2),  # output FASTQ pair
            ),
        )
    return input_output_pairs


@stage
class TrimAlignRNA(SequencingGroupStage):
    """
    Trim and align RNA-seq FASTQ reads with fastp and STAR
    """

    def expected_tmp_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Expect a pair of BAM and BAI files, one per set of input FASTQ files
        """
        return {
            suffix: sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('bam', 'bam'),
                ('bai', 'bam.bai'),
            ]
        }

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Expect a pair of CRAM and CRAI files, one per set of input FASTQ files
        """
        expected_outs = {
            suffix: sequencing_group.dataset.prefix() / 'cram' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('cram', 'cram'),
                ('crai', 'cram.crai'),
            ]
        }
        # Also include the temporary BAM and BAI files, but only if the CRAM and CRAI files don't exist
        if not (expected_outs['cram'].exists() and expected_outs['crai'].exists()):
            expected_outs.update(self.expected_tmp_outputs(sequencing_group))
        return expected_outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to align the input FASTQ files to the genome using STAR
        """
        jobs = []

        # Run trim
        input_fq_pairs = get_trim_inputs(sequencing_group)
        if not input_fq_pairs:
            return self.make_outputs(target=sequencing_group, error_msg='No FASTQ input found')
        assert isinstance(input_fq_pairs, FastqPairs)
        trimmed_fastq_pairs = []
        for fq_pair in input_fq_pairs:
            j, out_fqs = trim.trim(
                b=get_batch(),
                sequencing_group=sequencing_group,
                input_fq_pair=fq_pair,
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            if j:
                assert isinstance(j, Job)
                jobs.append(j)
            if not out_fqs or not isinstance(out_fqs, FastqPair):
                raise Exception(f'Error trimming FASTQs for {sequencing_group}')
            trimmed_fastq_pairs.append(out_fqs)

        # Run alignment
        trimmed_fastq_pairs = FastqPairs(trimmed_fastq_pairs)
        aligned_bam_dict = self.expected_tmp_outputs(sequencing_group)
        aligned_bam = BamPath(
            path=aligned_bam_dict['bam'],
            index_path=aligned_bam_dict['bai'],
        )
        aligned_cram_dict = self.expected_outputs(sequencing_group)
        aligned_cram = CramPath(
            path=aligned_cram_dict['cram'],
            index_path=aligned_cram_dict['crai'],
        )
        try:
            align_jobs = align_rna.align(
                b=get_batch(),
                fastq_pairs=trimmed_fastq_pairs,
                sample_name=sequencing_group.id,
                genome_prefix=get_config()['references']['star'].get('ref_dir'),
                mark_duplicates=True,
                output_bam=aligned_bam,
                output_cram=aligned_cram,
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            if align_jobs:
                assert isinstance(align_jobs, list)
                assert all([isinstance(j, Job) for j in align_jobs])
                jobs.extend(align_jobs)
        except Exception as e:
            logging.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
            raise Exception(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')

        # Create outputs and return jobs
        return self.make_outputs(sequencing_group, data=aligned_cram_dict, jobs=jobs)
