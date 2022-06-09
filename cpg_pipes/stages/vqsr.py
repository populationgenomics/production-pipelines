"""
Stage that performs AS-VQSR.
"""

import logging

from cpg_utils.config import get_config

from .. import Path
from ..pipeline import stage, CohortStage, StageInput, StageOutput
from ..jobs import vqsr
from ..targets import Cohort
from .joint_genotyping import JointGenotyping

logger = logging.getLogger(__file__)


@stage(required_stages=JointGenotyping)
class Vqsr(CohortStage):
    """
    Variant filtering of joint-called VCF.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expects to generate one site-only VCF.
        """
        h = cohort.alignment_inputs_hash()
        return self.tmp_prefix / f'{h}-siteonly.vcf.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        siteonly_vcf_path = inputs.as_path(
            stage=JointGenotyping, target=cohort, id='siteonly'
        )
        jobs = vqsr.make_vqsr_jobs(
            b=self.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            tmp_prefix=self.tmp_prefix,
            gvcf_count=len(cohort.get_samples()),
            output_vcf_path=self.expected_outputs(cohort),
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            overwrite=not get_config()['workflow'].get('self.check_intermediates'),
            scatter_count=get_config()['workflow'].get(
                'jc_intervals_num', vqsr.DEFAULT_INTERVALS_NUM
            ),
            sequencing_type=cohort.sequencing_type,
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
