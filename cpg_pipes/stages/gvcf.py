"""
Stage that generates a GVCF file.
"""

import logging

from .. import Path
from ..jobs import split_intervals, haplotype_caller
from ..pipeline.targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from .cram import CramStage

logger = logging.getLogger(__file__)


@stage(required_stages=CramStage, analysis_type='gvcf')
class GvcfStage(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples
    """
    hc_intervals = None

    def expected_result(self, sample: Sample) -> Path:
        """
        Generate a GVCF and corresponding TBI index
        """
        return sample.get_gvcf_path().path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Use function from the jobs module
        """
        hc_shards_num = self.pipeline_config.get('hc_shards_num', 1)
        if GvcfStage.hc_intervals is None and hc_shards_num > 1:
            GvcfStage.hc_intervals = split_intervals.get_intervals(
                b=self.b,
                refs=self.refs,
                scatter_count=hc_shards_num,
            )
        jobs = haplotype_caller.produce_gvcf(
            b=self.b,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            job_attrs=sample.get_job_attrs(),
            cram_path=sample.get_cram_path(),
            intervals=GvcfStage.hc_intervals,
            number_of_intervals=hc_shards_num,
            refs=self.refs,
            tmp_bucket=sample.dataset.get_tmp_bucket(),
            overwrite=not self.check_intermediates,
        )
        return self.make_outputs(
            sample,
            data=self.expected_result(sample), 
            jobs=jobs
        )
