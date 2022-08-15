"""
Stage that performs AS-VQSR.
"""

import logging

from cpg_utils import to_path
from cpg_utils.config import get_config

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

    def expected_outputs(self, cohort: Cohort):
        """
        Generate a site-only VCF.
        """
        h = cohort.alignment_inputs_hash()
        prefix = str(cohort.analysis_dataset.tmp_prefix() / self.name / h)
        return {
            'prefix': prefix,
            'siteonly': to_path(f'{prefix}-siteonly.vcf.gz'),
        }

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
            gvcf_count=len(cohort.get_samples()),
            out_path=self.expected_outputs(cohort)['siteonly'],
            tmp_prefix=to_path(self.expected_outputs(cohort)['prefix']),
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            scatter_count=get_config()['workflow'].get(
                'vqsr_intervals_num', vqsr.DEFAULT_INTERVALS_NUM
            ),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
