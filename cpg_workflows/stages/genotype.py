"""
Stage that generates a GVCF file.
"""
import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.workflows.filetypes import GvcfPath
from cpg_utils.workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)

from cpg_workflows.jobs import genotype
from cpg_workflows.jobs.happy import happy
from cpg_workflows.jobs.picard import vcf_qc

from .align import Align


@stage(required_stages=Align, analysis_type='gvcf')
class Genotype(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples (i.e. CRAM -> GVCF).
    """

    def expected_outputs(self, sample: Sample) -> dict:
        """
        Generate a GVCF and corresponding TBI index, as well as QC.
        """
        outs = {
            'gvcf': sample.make_gvcf_path().path,
        }
        if get_config()['workflow'].get('skip_qc', False) is False:
            qc_prefix = sample.dataset.prefix() / 'qc' / sample.id
            outs |= {
                'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
                'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
            }
        return outs

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        jobs = []
        gvcf_path = self.expected_outputs(sample)['gvcf']
        gvcf_jobs = genotype.genotype(
            b=self.b,
            output_path=gvcf_path,
            sample_name=sample.id,
            cram_path=sample.make_cram_path(),
            tmp_prefix=self.tmp_prefix / sample.id,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.extend(gvcf_jobs)
        if get_config()['workflow'].get('skip_qc', False) is False:
            qc_j = vcf_qc(
                b=self.b,
                vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(self.b),
                is_gvcf=True,
                job_attrs=self.get_job_attrs(sample),
                output_summary_path=self.expected_outputs(sample)['qc_summary'],
                output_detail_path=self.expected_outputs(sample)['qc_detail'],
                overwrite=not get_config()['workflow'].get('check_intermediates'),
            )
            if qc_j:
                if gvcf_jobs:
                    qc_j.depends_on(*gvcf_jobs)
                jobs.append(qc_j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)


@stage(required_stages=Genotype)
class GvcfHappy(SampleStage):
    """
    Run Happy to validate a GVCF for samples where a truth callset is available.
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
        return (
            sample.dataset.prefix()
            / 'qc'
            / 'gvcf'
            / 'hap.py'
            / f'{sample.id}.summary.csv'
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        gvcf_path = sample.make_gvcf_path()
        if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logging.warning(f'No GVCF found, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No GVCF found')

        jobs = happy(
            b=self.b,
            sample=sample,
            vcf_or_gvcf=sample.make_gvcf_path().resource_group(self.b),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sample),
            output_path=self.expected_outputs(sample),
        )

        if not jobs:
            return self.make_outputs(sample)
        else:
            return self.make_outputs(sample, self.expected_outputs(sample), jobs)
