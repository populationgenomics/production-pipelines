"""
Stage that generates a GVCF file.
"""

import logging
import hailtop.batch as hb
from cpg_utils.config import get_config

from .. import Path
from ..jobs import split_intervals, haplotype_caller
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

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Generate a GVCF and corresponding TBI index
        """
        return sample.get_gvcf_path().path

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
            intervals_j, intervals = split_intervals.get_intervals(
                b=self.b,
                intervals_path=get_config()['workflow'].get('intervals_path'),
                sequencing_type=self.cohort.sequencing_type,
                scatter_count=scatter_count,
                job_attrs=self.get_job_attrs(),
                output_prefix=self.tmp_prefix / 'intervals',
            )
            jobs.append(intervals_j)
            GenotypeSample.hc_intervals = intervals
        jobs.extend(
            haplotype_caller.produce_gvcf(
                b=self.b,
                output_path=self.expected_outputs(sample),
                sample_name=sample.id,
                cram_path=sample.get_cram_path(),
                intervals=GenotypeSample.hc_intervals,
                tmp_prefix=self.tmp_prefix / sample.id,
                overwrite=not get_config()['workflow'].get('self.check_intermediates'),
                job_attrs=self.get_job_attrs(sample),
            )
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
