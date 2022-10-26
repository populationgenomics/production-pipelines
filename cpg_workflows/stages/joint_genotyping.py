"""
Stage that performs joint genotyping of GVCFs using GATK.
"""
import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.workflows.filetypes import GvcfPath
from cpg_utils.workflows.workflow import (
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
            b=self.b,
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

        for job in jobs:
            assert job
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=JointGenotyping)
class JointVcfHappy(SampleStage):
    """
    Run Happy to validate validation samples in joint VCF
    """

    def expected_outputs(self, sample: Sample) -> Path | None:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        if sample.participant_id not in get_config().get('validation', {}).get(
            'sample_map', {}
        ):
            return None

        h = self.cohort.alignment_inputs_hash()
        return (
            self.cohort.analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / 'hap.py'
            / f'{h}-{sample.id}.summary.csv'
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        vcf_path = inputs.as_path(target=self.cohort, stage=JointGenotyping, key='vcf')

        jobs = happy(
            b=self.b,
            sample=sample,
            vcf_or_gvcf=self.b.read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
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
