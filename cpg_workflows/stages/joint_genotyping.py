"""
Stage that performs joint genotyping of GVCFs using GATK.
"""
import logging

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs import joint_genotyping
from cpg_workflows.workflow import (
    Cohort,
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    WorkflowError,
)
from .genotype import Genotype
from .. import get_batch
from ..resources import joint_calling_scatter_count


@stage(required_stages=Genotype)
class JointGenotyping(CohortStage):
    """
    Joint-calling of GVCFs together.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Generate a pVCF and a site-only VCF.
        """
        return {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.prefix / 'tmp'),
            'vcf': to_path(self.prefix / 'full.vcf.gz'),
            'siteonly': to_path(self.prefix / 'siteonly.vcf.gz'),
            'siteonly_part_pattern': str(
                self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        gvcf_by_sid = {
            sample.id: GvcfPath(
                inputs.as_path(target=sample, stage=Genotype, key='gvcf')
            )
            for sample in cohort.get_samples()
        }

        not_found_gvcfs: list[str] = []
        for sid, gvcf_path in gvcf_by_sid.items():
            if gvcf_path is None:
                logging.error(f'Joint genotyping: could not find GVCF for {sid}')
                not_found_gvcfs.append(sid)
        if not_found_gvcfs:
            raise WorkflowError(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting'
            )

        jobs = []
        vcf_path = self.expected_outputs(cohort)['vcf']
        siteonly_vcf_path = self.expected_outputs(cohort)['siteonly']
        scatter_count = joint_calling_scatter_count(len(cohort.get_samples()))
        out_siteonly_vcf_part_paths = [
            to_path(
                self.expected_outputs(cohort)['siteonly_part_pattern'].format(idx=idx)
            )
            for idx in range(scatter_count)
        ]

        jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
            b=get_batch(),
            out_vcf_path=vcf_path,
            out_siteonly_vcf_path=siteonly_vcf_path,
            tmp_bucket=to_path(self.expected_outputs(cohort)['tmp_prefix']),
            gvcf_by_sid=gvcf_by_sid,
            tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
            if get_config()['workflow'].get('use_gnarly', False)
            else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
            out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
        )
        jobs.extend(jc_jobs)
        for job in jobs:
            assert job
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
