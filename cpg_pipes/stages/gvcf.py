"""
Stage that generates a GVCF file.
"""

import logging
import hailtop.batch as hb

from .. import Path
from ..jobs import split_intervals, haplotype_caller
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from .cram import CramStage


logger = logging.getLogger(__file__)


@stage(required_stages=CramStage, analysis_type='gvcf')
class GvcfStage(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples
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
        hc_intervals_num = self.pipeline_config.get('hc_intervals_num', 1)
        jobs = []
        if GvcfStage.hc_intervals is None and hc_intervals_num > 1:
            intervals_j, intervals = split_intervals.get_intervals(
                b=self.b,
                refs=self.refs,
                sequencing_type=sample.sequencing_type,
                scatter_count=hc_intervals_num,
            )
            jobs.append(intervals_j)
            GvcfStage.hc_intervals = intervals
        jobs.extend(haplotype_caller.produce_gvcf(
            b=self.b,
            output_path=self.expected_outputs(sample),
            sample_name=sample.id,
            sequencing_type=sample.sequencing_type,
            job_attrs=sample.get_job_attrs(),
            cram_path=sample.get_cram_path(),
            intervals=GvcfStage.hc_intervals,
            number_of_intervals=hc_intervals_num,
            refs=self.refs,
            tmp_bucket=sample.dataset.get_tmp_bucket(),
            overwrite=not self.check_intermediates,
        ))
        return self.make_outputs(
            sample,
            data=self.expected_outputs(sample), 
            jobs=jobs
        )
