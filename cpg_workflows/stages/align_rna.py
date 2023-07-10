"""
Align RNA-seq reads to the genome using STAR.
"""

import logging
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import ExpectedResultT
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
)
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
    BamPath,
)
from cpg_workflows.stages.trim import Trim
from cpg_workflows.jobs import align_rna
import re


def get_alignment_inputs(sequencing_group: SequencingGroup, inputs: StageInput) -> FastqPairs:
    """
    Get the input FASTQ file pairs for alignment
    """
    input_fastqs = {
        k: v
        for k, v in inputs.as_dict(sequencing_group, Trim).items()
        if k.startswith('fastq_')
    }
    input_fastq_pairs_dict = {}
    for fq in input_fastqs:
        prefix = re.sub(r'_R[12]$', '', fq)
        if prefix not in input_fastq_pairs_dict:
            input_fastq_pairs_dict[prefix] = {}
        if fq.endswith('_R1'):
            input_fastq_pairs_dict[prefix]['R1'] = input_fastqs[fq]
        elif fq.endswith('_R2'):
            input_fastq_pairs_dict[prefix]['R2'] = input_fastqs[fq]
    input_fastq_pairs = FastqPairs([
        FastqPair(
            input_fastq_pairs_dict[prefix]['R1'],
            input_fastq_pairs_dict[prefix]['R2'],
        )
        for prefix in input_fastq_pairs
    ])
    return input_fastq_pairs


@stage(
    required_stages=Trim,
)
class AlignRNA(SequencingGroupStage):
    """
    Align RNA-seq FASTQ reads with STAR
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> ExpectedResultT:
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
        input_fastq_pairs = get_alignment_inputs(sequencing_group, inputs)
        _exp_out = self.expected_outputs(sequencing_group)
        output_bam = BamPath(
            path=_exp_out['bam'],
            index_path=_exp_out['bai'],
        )
        try:
            jobs = align_rna.align(
                b=get_batch(),
                input_fastq_pairs=input_fastq_pairs,
                output_bam=output_bam,
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