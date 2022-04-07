"""
VEP stage.
"""

from cpg_pipes import Path, utils
from cpg_pipes.jobs import vep
from cpg_pipes.pipeline import (
    stage,
    StageInput,
    StageOutput,
    CohortStage
)
from cpg_pipes.refdata import RefData
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
        samples_hash = utils.hash_sample_ids(cohort.get_sample_ids())
        return cohort.analysis_dataset.get_tmp_bucket() / 'vep' / f'{samples_hash}.ht'

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
            overwrite=not self.check_intermediates,
            scatter_count=self.pipeline_config.get(
                'vep_intervals_num', RefData.number_of_vep_intervals
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
