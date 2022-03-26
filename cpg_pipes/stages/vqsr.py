"""
Stage that performs AS-VQSR.
"""

import logging

from .. import Path
from .. import utils
from ..pipeline import stage, CohortStage, StageInput, StageOutput, Cohort
from ..jobs.vqsr import make_vqsr_jobs
from .joint_genotyping import JointGenotypingStage

logger = logging.getLogger(__file__)


@stage(required_stages=JointGenotypingStage)
class VqsrStage(CohortStage):
    """
    Variant filtering of joint-called VCF
    """
    def expected_result(self, cohort: Cohort) -> Path:
        """
        Expects to generate one site-only VCF
        """
        samples_hash = utils.hash_sample_ids(cohort.get_sample_ids())
        return (
            cohort.analysis_dataset.get_tmp_bucket() / 
            'vqsr' / 
            f'{samples_hash}-siteonly.vcf.gz'
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Use function defined in jobs module
        """
        siteonly_vcf_path = inputs.as_path(
            stage=JointGenotypingStage, target=cohort, id='siteonly'
        )

        tmp_vqsr_bucket = cohort.analysis_dataset.get_tmp_bucket() / 'vqsr'
        logger.info(f'Queueing VQSR job')
        expected_path = self.expected_result(cohort)
        jobs = make_vqsr_jobs(
            b=self.b,
            refs=self.refs,
            sequencing_type=cohort.get_sequencing_type(),
            input_vcf_or_mt_path=siteonly_vcf_path,
            work_bucket=tmp_vqsr_bucket,
            gvcf_count=len(cohort.get_samples()),
            output_vcf_path=expected_path,
            use_as_annotations=self.pipeline_config.get('use_as_vqsr', True),
            overwrite=not self.check_intermediates,
        )
        return self.make_outputs(cohort, data=expected_path, jobs=jobs)
