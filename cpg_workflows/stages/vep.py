"""
VEP stage.
"""

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    Cohort,
)

from cpg_workflows.jobs import vep
from .vqsr import Vqsr
from .. import get_batch


@stage(required_stages=[Vqsr])
class Vep(CohortStage):
    """
    Run VEP on a VCF.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a hail table.
        """
        h = cohort.alignment_inputs_hash()
        prefix = str(cohort.analysis_dataset.tmp_prefix() / self.name / h)
        return {
            'prefix': prefix,
            'ht': to_path(f'{prefix}.ht'),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        scatter_count = 50
        if len(cohort.get_samples()) > 300:
            scatter_count = 100
        if len(cohort.get_samples()) > 1000:
            scatter_count = 200

        jobs = vep.add_vep_jobs(
            get_batch(),
            vcf_path=inputs.as_path(cohort, stage=Vqsr, key='siteonly'),
            out_path=self.expected_outputs(cohort)['ht'],
            tmp_prefix=to_path(self.expected_outputs(cohort)['prefix']),
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
