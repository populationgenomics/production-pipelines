"""
VEP stage.
"""

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import vep
from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage

from ..resources import joint_calling_scatter_count
from .joint_genotyping import JointGenotyping
from .vqsr import Vqsr


@stage(required_stages=[Vqsr, JointGenotyping])
class Vep(MultiCohortStage):
    """
    Run VEP on a VCF.
    """

    def expected_outputs(self, multicohort: MultiCohort):
        """
        Expected to write a hail table.
        """
        return {'ht': self.prefix / 'vep.ht'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))
        input_siteonly_vcf_part_paths = [
            to_path(
                inputs.as_str(
                    stage=JointGenotyping,
                    target=multicohort,
                    key='siteonly_part_pattern',
                ).format(idx=idx),
            )
            for idx in range(scatter_count)
        ]

        jobs = vep.add_vep_jobs(
            get_batch(),
            input_vcfs=input_siteonly_vcf_part_paths,
            out_path=outputs['ht'],
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(multicohort, outputs, jobs)
