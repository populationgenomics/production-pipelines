"""
Stage that performs AS-VQSR.
"""

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import vqsr
from cpg_workflows.resources import joint_calling_scatter_count
from cpg_workflows.stages.joint_genotyping import JointGenotyping
from cpg_workflows.workflow import Cohort, CohortStage, StageInput, StageOutput, stage


@stage(required_stages=JointGenotyping)
class Vqsr(CohortStage):
    """
    Variant filtering of joint-called VCF.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Generate a site-only VCF.
        """
        return {'siteonly': self.prefix / 'siteonly.vcf.gz'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        This job takes the site-only VCF from the JointGenotyping stage with no prior
        multi-allelic variant splitting required
        VQSR itself splits out all variants
        """
        siteonly_vcf_path = inputs.as_path(stage=JointGenotyping, target=cohort, key='siteonly')
        outputs = self.expected_outputs(cohort)
        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=siteonly_vcf_path,
            scatter_count=joint_calling_scatter_count(len(cohort.get_sequencing_groups())),
            gvcf_count=len(cohort.get_sequencing_groups()),
            out_path=outputs['siteonly'],
            tmp_prefix=self.prefix / 'tmp',
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)
