"""
Stage that generates a GVCF file.
"""

import logging
import hailtop.batch as hb
from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from ..filetypes import GvcfPath
from ..jobs import haplotype_caller
from ..jobs.happy import happy
from ..jobs.picard import vcf_qc, get_intervals
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from .align import Align


logger = logging.getLogger(__file__)


@stage(required_stages=Align, analysis_type='gvcf')
class GenotypeSample(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples (i.e. CRAM -> GVCF).
    """

    hc_intervals: list[hb.Resource] | None = None

    def expected_outputs(self, sample: Sample) -> dict:
        """
        Generate a GVCF and corresponding TBI index, as well as QC.
        """
        qc_prefix = sample.dataset.prefix() / 'qc' / sample.id
        return {
            'gvcf': sample.make_gvcf_path().path,
            'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
            'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        scatter_count = get_config()['workflow'].get(
            'hc_intervals_num',
            haplotype_caller.DEFAULT_INTERVALS_NUM,
        )
        jobs = []
        if GenotypeSample.hc_intervals is None and scatter_count > 1:
            intervals_j, intervals = get_intervals(
                b=self.b,
                intervals_path=get_config()['workflow'].get('intervals_path'),
                scatter_count=scatter_count,
                job_attrs=self.get_job_attrs(),
                output_prefix=self.tmp_prefix / 'intervals',
            )
            jobs.append(intervals_j)
            GenotypeSample.hc_intervals = intervals
        gvcf_path = self.expected_outputs(sample)['gvcf']
        gvcf_jobs = haplotype_caller.produce_gvcf(
            b=self.b,
            output_path=gvcf_path,
            sample_name=sample.id,
            cram_path=sample.make_cram_path(),
            intervals=GenotypeSample.hc_intervals,
            tmp_prefix=self.tmp_prefix / sample.id,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.extend(gvcf_jobs)
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
        gvcf_path = sample.make_gvcf_path()
        if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(f'No GVCF found, skipping sample {sample}')
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
