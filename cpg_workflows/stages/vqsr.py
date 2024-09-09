"""
Stage that performs AS-VQSR.
"""

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import vqsr
from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage

from ..resources import joint_calling_scatter_count
from .joint_genotyping import JointGenotyping


@stage(required_stages=JointGenotyping)
class Vqsr(MultiCohortStage):
    """
    Variant filtering of joint-called VCF.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Generate a site-only VCF.
        """
        return {'siteonly': self.prefix / 'siteonly.vcf.gz'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        siteonly_vcf_path = inputs.as_path(stage=JointGenotyping, target=multicohort, key='siteonly')
        outputs = self.expected_outputs(multicohort)
        number_of_sgids = len(multicohort.get_sequencing_groups())

        intervals_path = None
        if config_retrieve(['workflow', 'intervals_path'], default=None):
            intervals_path = to_path(config_retrieve(['workflow', 'intervals_path']))

        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=siteonly_vcf_path,
            scatter_count=joint_calling_scatter_count(number_of_sgids),
            gvcf_count=number_of_sgids,
            out_path=outputs['siteonly'],
            tmp_prefix=self.tmp_prefix / 'tmp',
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=intervals_path,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)
