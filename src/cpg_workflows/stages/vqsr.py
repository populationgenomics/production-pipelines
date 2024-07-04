"""
Stage that performs AS-VQSR.
"""

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.jobs import vqsr
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)

from .. import get_batch
from ..resources import joint_calling_scatter_count
from .joint_genotyping import JointGenotyping


@stage(required_stages=JointGenotyping)
class Vqsr(CohortStage):
    """
    Variant filtering of joint-called VCF.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Generate a site-only VCF.
        """
        return {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.prefix / 'tmp'),
            'siteonly': self.prefix / 'siteonly.vcf.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        siteonly_vcf_path = inputs.as_path(stage=JointGenotyping, target=cohort, key='siteonly')

        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=siteonly_vcf_path,
            scatter_count=joint_calling_scatter_count(len(cohort.get_sequencing_groups())),
            gvcf_count=len(cohort.get_sequencing_groups()),
            out_path=self.expected_outputs(cohort)['siteonly'],
            tmp_prefix=to_path(self.expected_outputs(cohort)['tmp_prefix']),
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
