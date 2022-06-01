"""
Stage that performs joint genotyping of GVCFs using GATK.
"""

import logging

from .. import Path
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

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate a pVCF and a site-only VCF. Returns 2 outputs, thus not checking
        the SMDB, because the Analysis entry supports only single output.
        """
        h = cohort.alignment_inputs_hash()
        return {
            'vcf': self.tmp_bucket / f'{h}.vcf.gz',
            'siteonly': self.tmp_bucket / f'{h}-siteonly.vcf.gz',
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
            tmp_bucket=self.tmp_bucket,
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not self.check_intermediates,
            tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
            if self.pipeline_config.get('use_gnarly', False)
            else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
            scatter_count=self.pipeline_config.get(
                'jc_intervals_num', joint_genotyping.DEFAULT_INTERVALS_NUM
            ),
            sequencing_type=cohort.get_sequencing_type(),
            intervals_path=self.pipeline_config.get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
