"""
Stages that perform alignment QC on CRAM files.
"""

import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from cpg_pipes.jobs.happy import happy
from cpg_pipes.jobs.picard import vcf_qc
from cpg_pipes.pipeline import stage, SampleStage, StageInput, StageOutput, CohortStage
from cpg_pipes.targets import Sample, Cohort
from .genotype_sample import GenotypeSample
from .joint_genotyping import JointGenotyping

logger = logging.getLogger(__file__)


@stage
class GvcfQc(SampleStage):
    """
    GCVF QC using picard CollectVariantCallingMetrics
    Based on https://github.com/broadinstitute/warp/blob/d7e9c7de683dc6c70e6c32d580ba3a1898955f30/tasks/broad/Qc.wdl#L649
    """

    def expected_outputs(self, sample: Sample) -> dict:
        """
        One file (variant_calling_detail_metrics) is parsed by MultiQC
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L534-L538
        """
        prefix = sample.dataset.prefix() / 'qc' / sample.id
        return {
            'summary': to_path(f'{prefix}.variant_calling_summary_metrics'),
            'detail': to_path(f'{prefix}.variant_calling_detail_metrics'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        gvcf_path = sample.get_gvcf_path()
        if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(f'No GVCF found, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No GVCF found')

        j = vcf_qc(
            b=self.b,
            vcf_or_gvcf=gvcf_path.resource_group(self.b),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sample),
            output_summary_path=self.expected_outputs(sample)['summary'],
            output_detail_path=self.expected_outputs(sample)['detail'],
        )

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[j])


@stage
class JointVcfQc(CohortStage):
    """
    GCVF QC using picard CollectVariantCallingMetrics
    Based on https://github.com/broadinstitute/warp/blob/d7e9c7de683dc6c70e6c32d580ba3a1898955f30/tasks/broad/Qc.wdl#L649
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        One file (variant_calling_detail_metrics) is parsed by MultiQC
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L534-L538
        """
        h = self.cohort.alignment_inputs_hash()
        prefix = self.cohort.analysis_dataset.prefix() / 'qc' / 'jc' / 'picard'
        return {
            'summary': to_path(f'{prefix}.variant_calling_summary_metrics'),
            'detail': to_path(f'{prefix}.variant_calling_detail_metrics'),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, id='vcf')

        j = vcf_qc(
            b=self.b,
            vcf_or_gvcf=self.b.read_input_group(
                **{
                    'vcf': str(vcf_path),
                    'vcf.tbi': str(vcf_path) + '.tbi',
                }
            ),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(cohort),
            output_summary_path=self.expected_outputs(cohort)['summary'],
            output_detail_path=self.expected_outputs(cohort)['detail'],
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage
class GvcfHappy(SampleStage):
    """
    Run Happy to validate GCVF of validation samples
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        return sample.dataset.prefix() / 'qc' / f'{sample.id}.summary.csv'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        gvcf_path = sample.get_gvcf_path()
        if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(f'No GVCF found, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No GVCF found')

        jobs = happy(
            b=self.b,
            sample=sample,
            vcf_or_gvcf=sample.get_gvcf_path().resource_group(self.b),
            is_gvcf=True,
            seqtype=self.cohort.sequencing_type,
            job_attrs=self.get_job_attrs(sample),
            output_path=self.expected_outputs(sample),
        )

        if not jobs:
            return self.make_outputs(sample)
        else:
            return self.make_outputs(sample, self.expected_outputs(sample), jobs)


@stage
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
            seqtype=self.cohort.sequencing_type,
            job_attrs=self.get_job_attrs(sample),
            output_path=self.expected_outputs(sample),
        )
        if not jobs:
            return self.make_outputs(sample)
        else:
            return self.make_outputs(sample, self.expected_outputs(sample), jobs)
