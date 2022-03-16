"""
Stage that generates a GVCF file.
"""

import logging

from cpg_pipes.jobs import split_intervals, haplotype_caller
from cpg_pipes.pipeline.analysis import AnalysisType
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput
from cpg_pipes.stages.cram import CramStage

logger = logging.getLogger(__file__)


@stage(required_stages=CramStage, sm_analysis_type=AnalysisType.GVCF)
class GvcfStage(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples
    """
    hc_intervals = None

    def expected_result(self, sample: Sample):
        """
        Generate a GVCF and corresponding TBI index
        """
        return sample.get_gvcf_path().path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Use function from the jobs module
        """
        hc_shards_num = self.pipe.config.get('hc_shards_num', 1)
        if GvcfStage.hc_intervals is None and hc_shards_num > 1:
            GvcfStage.hc_intervals = split_intervals.get_intervals(
                b=self.pipe.b,
                scatter_count=hc_shards_num,
            )
        gvcf_job = haplotype_caller.produce_gvcf(
            b=self.pipe.b,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
            cram_path=sample.get_cram_path(),
            intervals=GvcfStage.hc_intervals,
            number_of_intervals=hc_shards_num,
            tmp_bucket=self.pipe.tmp_bucket,
            overwrite=not self.pipe.check_intermediates,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
        )
        return self.make_outputs(
            sample,
            data=self.expected_result(sample), 
            jobs=[gvcf_job]
        )
