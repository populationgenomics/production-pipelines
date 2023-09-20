"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""

import logging
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
)
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.stages.trim_align import TrimAlignRNA
from cpg_workflows.jobs import count
from cpg_workflows.jobs import bam_to_cram


@stage(
    required_stages=TrimAlignRNA,
)
class Count(SequencingGroupStage):
    """
    Count reads with featureCounts.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate a text file output containing read counts.
        """
        return {
            'count': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count',
            'summary': sequencing_group.dataset.prefix() / 'count' / f'{sequencing_group.id}.count.summary',
        }
    
    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to count the reads with featureCounts.
        """
        cram_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'cram')
        crai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'crai')
        potential_bam_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam'
        potential_bai_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam.bai'
        input_bam = BamPath(potential_bam_path, potential_bai_path)
        output_path = self.expected_outputs(sequencing_group)['count']
        summary_path = self.expected_outputs(sequencing_group)['summary']
        jobs: list[Job] = []
        if (
            can_reuse(output_path, overwrite=sequencing_group.forced) and
            can_reuse(summary_path, overwrite=sequencing_group.forced)
        ):
            logging.info(f'Reusing existing count output {output_path}')
            return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
        elif potential_bam_path.exists() and potential_bai_path.exists():
            logging.info(f'Using existing BAM {potential_bam_path}')
        else:
            j = bam_to_cram.cram_to_bam(
                b=get_batch(),
                input_cram=CramPath(cram_path, crai_path),
                output_bam=input_bam,
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
            )
            jobs.append(j)
        jobs.append(
            count.count(
                b=get_batch(),
                input_bam=input_bam,
                output_path=output_path,
                summary_path=summary_path,
                sample_name=sequencing_group.id,
                job_attrs=self.get_job_attrs(sequencing_group),
            )
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
