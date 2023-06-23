"""
Stage that generates a GVCF file.
"""

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows.jobs import genotype
from cpg_workflows.workflow import (
    SequencingGroup,
    stage,
    StageInput,
    StageOutput,
    SequencingGroupStage,
)
from .align import Align
from .. import get_batch


@stage(
    required_stages=Align,
    analysis_type='gvcf',
    analysis_keys=['gvcf'],
)
class Genotype(SequencingGroupStage):
    """
    Use HaplotypeCaller to genotype individual sequencing groups (i.e. CRAM -> GVCF).
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate a GVCF and corresponding TBI index.
        """
        return {
            'gvcf': sequencing_group.make_gvcf_path().path,
            'tbi': sequencing_group.make_gvcf_path().tbi_path,
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        jobs = genotype.genotype(
            b=get_batch(),
            output_path=self.expected_outputs(sequencing_group)['gvcf'],
            sequencing_group_name=sequencing_group.id,
            cram_path=sequencing_group.make_cram_path(),
            tmp_prefix=self.tmp_prefix / sequencing_group.id,
            overwrite=sequencing_group.forced,
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        return self.make_outputs(
            sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs
        )
