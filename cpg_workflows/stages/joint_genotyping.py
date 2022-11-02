"""
Stage that performs joint genotyping of GVCFs using GATK.
"""
import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.workflow import (
    Sample,
    Cohort,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    CohortStage,
    WorkflowError,
)
from cpg_workflows.jobs.happy import happy
from cpg_workflows.jobs.picard import vcf_qc
from cpg_workflows.jobs import joint_genotyping
from .genotype import Genotype
from .. import get_batch


@stage(required_stages=Genotype)
class JointGenotyping(CohortStage):
    """
    Joint-calling of GVCFs together.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Generate a pVCF and a site-only VCF.
        """
        h = cohort.alignment_inputs_hash()
        prefix = str(cohort.analysis_dataset.prefix() / self.name / h)
        qc_prefix = cohort.analysis_dataset.prefix() / 'qc' / 'jc' / h / 'picard'
        return {
            'prefix': prefix,
            'vcf': to_path(f'{prefix}.vcf.gz'),
            'siteonly': to_path(f'{prefix}-siteonly.vcf.gz'),
            'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
            'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
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

        jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
            b=get_batch(),
            out_vcf_path=vcf_path,
            out_siteonly_vcf_path=siteonly_vcf_path,
            tmp_bucket=to_path(self.expected_outputs(cohort)['prefix']),
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
            if get_config()['workflow'].get('use_gnarly', False)
            else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        jobs.extend(jc_jobs)

        qc_j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=get_batch().read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
                }
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(cohort),
            output_summary_path=self.expected_outputs(cohort)['qc_summary'],
            output_detail_path=self.expected_outputs(cohort)['qc_detail'],
        )
        if qc_j:
            qc_j.depends_on(*jc_jobs)
            jobs.append(qc_j)

        for job in jobs:
            assert job
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
