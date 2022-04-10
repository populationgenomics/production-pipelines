"""
Stage that performs AS-VQSR.
"""

import logging

from .. import Path
from ..pipeline import stage, CohortStage, StageInput, StageOutput
from ..jobs.vqsr import make_vqsr_jobs
from ..targets import Cohort
from .joint_genotyping import JointGenotypingStage

logger = logging.getLogger(__file__)


@stage(required_stages=JointGenotypingStage)
class VqsrStage(CohortStage):
    """
    Variant filtering of joint-called VCF.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expects to generate one site-only VCF.
        """
        h = cohort.alignment_inputs_hash()
        return self.tmp_bucket / f'{h}-siteonly.vcf.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs.
        """
        siteonly_vcf_path = inputs.as_path(
            stage=JointGenotypingStage, target=cohort, id='siteonly'
        )
        jobs = make_vqsr_jobs(
            b=self.b,
            refs=self.refs,
            sequencing_type=cohort.get_sequencing_type(),
            input_vcf_or_mt_path=siteonly_vcf_path,
            tmp_bucket=self.tmp_bucket,
            gvcf_count=len(cohort.get_samples()),
            output_vcf_path=self.expected_outputs(cohort),
            use_as_annotations=self.pipeline_config.get('use_as_vqsr', True),
            overwrite=not self.check_intermediates,
            scatter_count=self.pipeline_config.get('jc_intervals_num'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
