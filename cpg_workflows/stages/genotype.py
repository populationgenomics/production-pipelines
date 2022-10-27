"""
Stage that generates a GVCF file.
"""
import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config

from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.workflow import (
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

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Generate a GVCF and corresponding TBI index.
        """
        return {
            'gvcf': sample.make_gvcf_path().path,
            'tbi': sample.make_gvcf_path().tbi_path,
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        jobs = genotype.genotype(
            b=self.b,
            output_path=self.expected_outputs(sample)['gvcf'],
            sample_name=sample.id,
            cram_path=sample.make_cram_path(),
            tmp_prefix=self.tmp_prefix / sample.id,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
