"""
Stage that performs joint genotyping of GVCFs using GATK.
"""

import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from cpg_pipes.jobs.happy import happy
from cpg_pipes.jobs.picard import vcf_qc
from cpg_pipes.filetypes import GvcfPath
from cpg_pipes.jobs import joint_genotyping
from cpg_pipes.targets import Cohort, Sample
from cpg_pipes.pipeline import (
    stage,
    CohortStage,
    StageInput,
    StageOutput,
    PipelineError,
    SampleStage,
)
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
        qc_prefix = self.cohort.analysis_dataset.prefix() / 'qc' / 'jc' / h / 'picard'
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
                inputs.as_path(target=sample, stage=GenotypeSample, id='gvcf')
            )
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

        jobs = []
        vcf_path = self.expected_outputs(cohort)['vcf']
        siteonly_vcf_path = self.expected_outputs(cohort)['siteonly']

        jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
            b=self.b,
            out_vcf_path=vcf_path,
            out_siteonly_vcf_path=siteonly_vcf_path,
            tmp_bucket=to_path(self.expected_outputs(cohort)['prefix']),
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
            if get_config()['workflow'].get('use_gnarly', False)
            else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
            scatter_count=get_config()['workflow'].get(
                'jc_intervals_num', joint_genotyping.DEFAULT_INTERVALS_NUM
            ),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        jobs.extend(jc_jobs)

        qc_j = vcf_qc(
            b=self.b,
            vcf_or_gvcf=self.b.read_input_group(
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

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=JointGenotyping)
class JointVcfHappy(SampleStage):
    """
    Run Happy to validate validation samples in joint VCF
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        h = self.cohort.alignment_inputs_hash()
        prefix = self.cohort.analysis_dataset.prefix() / 'qc' / 'jc' / 'happy'
        return prefix / f'{h}-{sample.id}.summary.csv'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        vcf_path = inputs.as_path(target=self.cohort, stage=JointGenotyping, id='vcf')

        jobs = happy(
            b=self.b,
            sample=sample,
            vcf_or_gvcf=self.b.read_input_group(
                **{
                    'vcf': str(vcf_path),
                    'vcf.tbi': str(vcf_path) + '.tbi',
                }
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(sample),
            output_path=self.expected_outputs(sample),
        )
        if not jobs:
            return self.make_outputs(sample)
        else:
            return self.make_outputs(sample, self.expected_outputs(sample), jobs)
