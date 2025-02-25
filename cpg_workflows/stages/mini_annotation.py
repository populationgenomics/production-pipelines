#!/usr/bin/env python3

"""
This file describes a workflow to apply minimal annotations to a callset
Starting from first principles -

- single-sample VCFs
- publicly available reference data
- no specialized tools (or just tools that can be built from source easily)
"""

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import SequencingGroup, MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage, DatasetStage, SequencingGroupStage
from cpg_workflows.jobs.bcftools import strip_gvcf_to_vcf


@stage
class StripSingleSampleGvcf(SequencingGroupStage):
    """
    Merge all the single-sample VCFs into a single VCF
    do this by finding the gVCFs and stripping them to VCF
    This isn't perfect, it's just to stress test this type of workflow
    We get the same results by condensing the combiner output to joint-VCF
    """

    def expected_outputs(self, sg: SequencingGroup) -> Path:
        return self.tmp_prefix / f'{sg.id}_stripped.vcf.bgz'

    def queue_jobs(self, sg: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        Take the gVCF for this sample and strip it to a VCF (kinda...)
        """

        outputs = self.expected_outputs(sg)

        job = strip_gvcf_to_vcf(
            gvcf=sg.gvcf,
            output=str(outputs),
        )
        return self.make_outputs(sg, data=outputs, jobs=job)


