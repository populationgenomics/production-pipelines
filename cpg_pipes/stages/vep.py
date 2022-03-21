"""
VEP stage.
"""

from cpg_pipes import Path
from cpg_pipes.jobs import vep
from cpg_pipes.pipeline import (
    Cohort,
    stage,
    StageInput,
    StageOutput,
    CohortStage
)
from cpg_pipes.stages.joint_genotyping import JointGenotypingStage


@stage(required_stages=[JointGenotypingStage])
class VepStage(CohortStage):
    """
    Run VEP on a VCF
    """
    def expected_result(self, cohort: Cohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return cohort.analysis_dataset.get_tmp_bucket() / 'vep' / 'cohort-vep.vcf.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Uses jobs vep() function.
        """
        j = vep.vep(
            self.b,
            vcf_path=inputs.as_path(cohort, stage=JointGenotypingStage, id='vcf'),
            refs=self.refs,
            out_vcf_path=self.expected_result(cohort),
        )
        return self.make_outputs(cohort, self.expected_result(cohort), [j])
