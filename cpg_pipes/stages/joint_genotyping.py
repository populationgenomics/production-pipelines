"""
Stage that performs joint genotyping of GVCFs using GATK.
"""

import logging

from cpg_pipes.storage import Path, to_path

from cpg_pipes import utils
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.pipeline.analysis import GvcfPath
from cpg_pipes.pipeline.dataset import Cohort
from cpg_pipes.pipeline.pipeline import PipelineError, stage, CohortStage
from cpg_pipes.pipeline.stage import StageInput, StageOutput
from cpg_pipes.stages.gvcf import GvcfStage

logger = logging.getLogger(__file__)


@stage(required_stages=GvcfStage)
class JointGenotypingStage(CohortStage):
    """
    Joint-calling of GVCFs together.
    """
    def expected_result(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate a pVCF and a site-only VCF. Returns 2 outputs, thus not checking
        the SMDB, because the Analysis entry supports only single output.
        """
        samples_hash = utils.hash_sample_ids(cohort.get_all_sample_ids())
        expected_jc_vcf_path = (
            cohort.analysis_dataset.get_tmp_bucket() / 
            'joint_calling' / 
            f'{samples_hash}.vcf.gz'
        )
        return {
            'vcf': expected_jc_vcf_path,
            'siteonly': to_path(
                str(expected_jc_vcf_path).replace('.vcf.gz', '-siteonly.vcf.gz')
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Use function defined in jobs module
        """
        gvcf_by_sid = {
            sample.id: 
                GvcfPath(inputs.as_path(target=sample, stage=GvcfStage)) 
            for sample in cohort.get_all_samples()
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

        jobs = make_joint_genotyping_jobs(
            b=self.b,
            out_vcf_path=self.expected_result(cohort)['vcf'],
            out_siteonly_vcf_path=self.expected_result(cohort)['siteonly'],
            samples=cohort.get_all_samples(),
            genomicsdb_bucket=cohort.analysis_dataset.get_bucket() / 'genomicsdbs',
            tmp_bucket=cohort.analysis_dataset.get_tmp_bucket(),
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not self.check_intermediates,
            depends_on=inputs.get_jobs(),
            smdb=self.smdb,
            tool=JointGenotyperTool.GnarlyGenotyper 
            if self.pipeline_config.get('use_gnarly', False) 
            else JointGenotyperTool.GenotypeGVCFs,
            dry_run=self.dry_run,
        )
        return self.make_outputs(
            cohort, 
            data=self.expected_result(cohort), 
            jobs=jobs
        )
