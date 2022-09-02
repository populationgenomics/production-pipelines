"""
Stage that generates a GVCF file.
"""

import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)

from jobs import genotype
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
        qc_prefix = sample.dataset.prefix() / 'qc' / sample.id
        return {
            'gvcf': sample.make_gvcf_path().path,
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        gvcf_path = self.expected_outputs(sample)['gvcf']
        jobs = genotype.genotype(
            b=self.b,
            output_path=gvcf_path,
            sample_name=sample.id,
            cram_path=sample.make_cram_path(),
            tmp_prefix=self.tmp_prefix / sample.id,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
