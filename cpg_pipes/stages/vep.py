"""
VEP stage.
"""

from cpg_pipes import Path
from cpg_pipes.jobs import vep
from cpg_pipes.pipeline import (
    stage,
    StageInput,
    StageOutput,
    CohortStage
)
from cpg_pipes.stages.vqsr import VqsrStage
from cpg_pipes.targets import Cohort


@stage(required_stages=[VqsrStage])
class VepStage(CohortStage):
    """
    Run VEP on a VCF
    """
    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return cohort.analysis_dataset.get_tmp_bucket() / 'vep' / 'vep.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Uses jobs vep() function.
        """
        jobs = vep.vep(
            self.b,
            vcf_path=inputs.as_path(cohort, stage=VqsrStage),
            refs=self.refs,
            sequencing_type=cohort.get_sequencing_type(),
            hail_billing_project=self.hail_billing_project,
            hail_bucket=self.hail_bucket,
            out_path=self.expected_outputs(cohort),
            tmp_bucket=self.tmp_bucket,
            scatter_count=100,
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
