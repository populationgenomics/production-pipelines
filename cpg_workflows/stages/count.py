"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""

import logging
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
)
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.stages.align_rna import AlignRNA
from cpg_workflows.jobs import count


@stage(
    required_stages=AlignRNA,
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
        bam_path = inputs.as_path(sequencing_group, AlignRNA, 'bam')
        bai_path = inputs.as_path(sequencing_group, AlignRNA, 'bai')
        j = count.count(
            b=get_batch(),
            input_bam=BamPath(bam_path, bai_path),
            output_path=self.expected_outputs(sequencing_group)['count'],
            summary_path=self.expected_outputs(sequencing_group)['summary'],
            sample_name=sequencing_group.id,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=j)