"""
VEP stage.
"""

from cpg_pipes import Path
from cpg_pipes.jobs import vep
from cpg_pipes.pipeline import stage, StageInput, StageOutput, CohortStage
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.targets import Cohort


@stage(required_stages=[Vqsr])
class Vep(CohortStage):
    """
    Run VEP on a VCF.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expected to write a hail table.
        """
        h = cohort.alignment_inputs_hash()
        return self.tmp_bucket / f'{h}.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs.
        """
        jobs = vep.vep(
            self.b,
            vcf_path=inputs.as_path(cohort, stage=Vqsr),
            refs=self.refs,
            hail_billing_project=self.hail_billing_project,
            hail_bucket=self.hail_bucket,
            out_path=self.expected_outputs(cohort),
            tmp_bucket=self.tmp_bucket,
            overwrite=not self.check_intermediates,
            scatter_count=self.pipeline_config.get('vep_intervals_num'),
            sequencing_type=cohort.get_sequencing_type(),
            intervals_path=self.pipeline_config.get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
