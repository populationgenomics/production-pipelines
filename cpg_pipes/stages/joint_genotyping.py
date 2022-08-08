"""
Stage that performs joint genotyping of GVCFs using GATK.
"""

import logging

from cpg_utils import to_path
from cpg_utils.config import get_config

from ..types import GvcfPath
from ..jobs import joint_genotyping
from ..targets import Cohort
from ..pipeline import stage, CohortStage, StageInput, StageOutput, PipelineError
from .genotype_sample import GenotypeSample

logger = logging.getLogger(__file__)


@stage(required_stages=GenotypeSample)
class JointGenotyping(CohortStage):
    """
    Joint-calling of GVCFs together.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Generate a pVCF and a site-only VCF.
        """
        h = cohort.alignment_inputs_hash()
        prefix = str(cohort.analysis_dataset.tmp_prefix() / self.name / h)
        return {
            'prefix': prefix,
            'vcf': to_path(f'{prefix}.vcf.gz'),
            'siteonly': to_path(f'{prefix}-siteonly.vcf.gz'),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        gvcf_by_sid = {
            sample.id: GvcfPath(inputs.as_path(target=sample, stage=GenotypeSample))
            for sample in cohort.get_samples()
        }

        not_found_gvcfs: list[str] = []
        for sid, gvcf_path in gvcf_by_sid.items():
            if gvcf_path is None:
                logger.error(f'Joint genotyping: could not find GVCF for {sid}')
                not_found_gvcfs.append(sid)
        if not_found_gvcfs:
            raise PipelineError(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting'
            )

        jobs = joint_genotyping.make_joint_genotyping_jobs(
            b=self.b,
            out_vcf_path=self.expected_outputs(cohort)['vcf'],
            out_siteonly_vcf_path=self.expected_outputs(cohort)['siteonly'],
            tmp_bucket=to_path(self.expected_outputs(cohort)['prefix']),
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
            if get_config()['workflow'].get('use_gnarly', False)
            else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
            scatter_count=get_config()['workflow'].get(
                'jc_intervals_num', joint_genotyping.DEFAULT_INTERVALS_NUM
            ),
            sequencing_type=self.cohort.sequencing_type,
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
