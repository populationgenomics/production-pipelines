"""
Stage that generates a GVCF file.
"""

import logging
import hailtop.batch as hb

from .. import Path
from ..jobs import split_intervals, haplotype_caller
from ..refdata import RefData
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

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Use function from the jobs module
        """
        scatter_count = self.pipeline_config.get(
            'hc_intervals_num',
            RefData.number_of_haplotype_caller_intervals,
        )
        jobs = []
        if GenotypeSample.hc_intervals is None and scatter_count > 1:
            intervals_j, intervals = split_intervals.get_intervals(
                b=self.b,
                refs=self.refs,
                intervals_path=self.pipeline_config.get('intervals_path'),
                sequencing_type=sample.sequencing_type,
                scatter_count=scatter_count,
                job_attrs=self.get_job_attrs(),
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
                refs=self.refs,
                tmp_bucket=self.tmp_bucket / sample.id,
                overwrite=not self.check_intermediates,
                job_attrs=self.get_job_attrs(sample),
            )
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
