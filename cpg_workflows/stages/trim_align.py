"""
Align RNA-seq reads to the genome using STAR.
"""

import logging
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
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
    CramPath
)
from cpg_workflows.jobs import trim
from cpg_workflows.jobs import align_rna
from cpg_workflows.jobs import markdups
from cpg_workflows.jobs import bam_to_cram
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
        input_r1_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r1)))
        input_r2_bn = re.sub('.f(ast)?q.gz', '', basename(str(pair.r2)))
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
        Expect a pair of CRAM and CRAI files, one per set of input FASTQ files
        """
        return {
            suffix: sequencing_group.dataset.prefix() / 'cram' / f'{sequencing_group.id}.{extension}'
            for suffix, extension in [
                ('cram', 'cram'),
                ('crai', 'cram.crai'),
            ]
        }
    
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

    
    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to align the input FASTQ files to the genome using STAR
        """
        jobs = []

        # Run trim
        input_fq_pairs = get_trim_inputs(sequencing_group)
        if not input_fq_pairs:
            if get_config()['workflow'].get('skip_sgs_with_missing_input'):
                logging.error(f'No FASTQ inputs, skipping sample {sequencing_group}')
                sequencing_group.active = False
                return self.make_outputs(sequencing_group, skipped=True)  # return empty output
            else:
                return self.make_outputs(
                    target=sequencing_group, error_msg=f'No FASTQ input found'
                )
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
                if out_fqs:
                    assert isinstance(out_fqs, FastqPair)
                    trimmed_fastq_pairs.append(out_fqs)
            elif out_fqs:
                # If the job was skipped due to existing trimmed FASTQs,
                # we need to localise them
                assert isinstance(out_fqs, FastqPair)
                ex_fq = get_batch().read_input_group(
                    r1=str(out_fqs.r1),
                    r2=str(out_fqs.r2),
                )
                trimmed_fastq_pairs.append(FastqPair(
                    r1=ex_fq['r1'],
                    r2=ex_fq['r2'],
                ))

        # Run alignment
        trimmed_fastq_pairs = FastqPairs(trimmed_fastq_pairs)
        aligned_bam = None
        try:
            align_jobs, bam_out = align_rna.align(
                b=get_batch(),
                fastq_pairs=trimmed_fastq_pairs,
                sample_name=sequencing_group.id,
                genome_prefix=get_config()['references'].get('star_ref_dir'),
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            if align_jobs:
                assert isinstance(align_jobs, list)
                assert all([isinstance(j, Job) for j in align_jobs])
                assert isinstance(bam_out, ResourceGroup)
                jobs.extend(align_jobs)
                aligned_bam = bam_out
            else:
                # If the job was skipped due to existing BAMs,
                # we need to localise them
                assert isinstance(bam_out, BamPath)
                aligned_bam = get_batch().read_input_group(
                    bam=str(bam_out.path),
                    bai=str(bam_out.index_path),
                )
        except Exception as e:
            logging.error(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')
            raise Exception(f'Error aligning RNA-seq reads for {sequencing_group}: {e}')

        # Run mark duplicates
        _exp_tmp_out = self.expected_tmp_outputs(sequencing_group)
        output_bam = BamPath(
            path=_exp_tmp_out['bam'],
            index_path=_exp_tmp_out['bai'],
        )
        j, mkdup_bam_out = markdups.markdup(
            b=get_batch(),
            input_bam=aligned_bam,
            output_bam=output_bam,
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=sequencing_group.forced,
        )
        if j:
            assert isinstance(j, Job)
            assert isinstance(mkdup_bam_out, ResourceGroup)
            jobs.append(j)
            mkdup_bam = mkdup_bam_out
        else:
            # If the job was skipped due to an existing BAM,
            # we need to localise it
            assert isinstance(mkdup_bam_out, BamPath)
            mkdup_bam = get_batch().read_input_group(
                bam=str(mkdup_bam_out.path),
                bai=str(mkdup_bam_out.index_path),
            )

        # Convert BAM to CRAM
        _exp_out = self.expected_outputs(sequencing_group)
        output_cram = CramPath(
            path=_exp_out['cram'],
            index_path=_exp_out['crai'],
        )
        j, mkdup_cram = bam_to_cram.bam_to_cram(
            b=get_batch(),
            input_bam=mkdup_bam,
            output_cram=output_cram,
            extra_label='mkdup',
            job_attrs=self.get_job_attrs(sequencing_group),
            overwrite=sequencing_group.forced,
        )
        if j:
            assert isinstance(j, Job)
            assert isinstance(mkdup_cram, ResourceGroup)
            jobs.append(j)
        else:
            assert isinstance(mkdup_cram, CramPath)

        # Create outputs and return jobs
        return self.make_outputs(
            sequencing_group,
            data=_exp_out,
            jobs=jobs
        )
