"""
VEP stage.
"""

from cpg_utils import to_path
from cpg_workflows.jobs import vep
from cpg_workflows.resources import joint_calling_scatter_count
from cpg_workflows.stages.joint_genotyping import JointGenotyping
from cpg_workflows.targets import Cohort
from cpg_workflows.workflow import CohortStage, StageInput, StageOutput, stage


@stage(required_stages=JointGenotyping)
class Vep(CohortStage):
    """
    Run VEP on a VCF.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a hail table.
        """
        return {'ht': self.prefix / 'vep.ht'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(cohort)
        scatter_count = joint_calling_scatter_count(len(cohort.get_sequencing_groups()))
        input_vcf_paths = [
            to_path(
                inputs.as_str(
                    stage=JointGenotyping,
                    target=cohort,
                    key='siteonly_split_part_pattern',
                ).format(idx=idx),
            )
            for idx in range(scatter_count)
        ]

        jobs = vep.add_vep_jobs(
            input_vcf_paths=input_vcf_paths,
            out_path=outputs['ht'],
            tmp_prefix=self.prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, outputs, jobs)
