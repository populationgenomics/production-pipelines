"""
Stages that perform alignment QC on CRAM files.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config

from cpg_pipes.jobs.happy import happy
from cpg_pipes.pipeline import stage, SampleStage, StageInput, StageOutput, CohortStage
from cpg_pipes.targets import Sample
from .joint_genotyping import JointGenotyping

logger = logging.getLogger(__file__)


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
