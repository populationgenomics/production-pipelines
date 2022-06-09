"""
VEP stage.
"""
from cpg_utils.config import get_config

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

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        jobs = vep.vep_jobs(
            self.b,
            vcf_path=inputs.as_path(cohort, stage=Vqsr),
            out_path=self.expected_outputs(cohort),
            tmp_bucket=self.tmp_bucket,
            overwrite=not get_config()['workflow'].get('self.check_intermediates'),
            scatter_count=get_config()['workflow'].get('vep_intervals_num', vep.DEFAULT_INTERVALS_NUM),
            sequencing_type=get_config()['workflow']['sequencing_type'],
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
