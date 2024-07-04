"""
VEP stage.
"""

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import vep
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)

from ..resources import joint_calling_scatter_count
from .joint_genotyping import JointGenotyping
from .vqsr import Vqsr


@stage(required_stages=[Vqsr, JointGenotyping])
class Vep(CohortStage):
    """
    Run VEP on a VCF.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a hail table.
        """
        return {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.prefix / 'tmp'),
            'ht': self.prefix / 'vep.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        scatter_count = joint_calling_scatter_count(len(cohort.get_sequencing_groups()))
        input_siteonly_vcf_part_paths = [
            to_path(
                inputs.as_str(
                    stage=JointGenotyping,
                    target=cohort,
                    key='siteonly_part_pattern',
                ).format(idx=idx),
            )
            for idx in range(scatter_count)
        ]

        jobs = vep.add_vep_jobs(
            get_batch(),
            input_siteonly_vcf_path=inputs.as_path(cohort, stage=Vqsr, key='siteonly'),
            input_siteonly_vcf_part_paths=input_siteonly_vcf_part_paths,
            out_path=self.expected_outputs(cohort)['ht'],
            tmp_prefix=to_path(self.expected_outputs(cohort)['tmp_prefix']),
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), jobs)
