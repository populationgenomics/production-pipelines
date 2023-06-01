#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""


from typing import Any

from cpg_utils import to_path, Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    SampleStage,
    Sample,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
    Cohort,
    Dataset,
    get_workflow,
)
from cpg_workflows.stages.seqr_loader import AnnotateDataset, _sg_vcf_meta
from cpg_workflows.jobs.validation import validation_mt_to_vcf_job
from .. import get_batch


@stage(
    required_stages=AnnotateDataset,
    analysis_type='custom',
    update_analysis_meta=_sg_vcf_meta,
    analysis_keys=['vcf'],
)
class ValidationMtToVcf(SampleStage):

    def expected_outputs(self, sample: Sample):

        return {
            'vcf': (
                sample.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sample.id}.vcf.bgz'
            ),
            'index': (
                sample.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sample.id}.vcf.bgz.tbi'
            )
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        # only run this on validation samples
        if sample.dataset.name != 'validation':
            return None

        # get the input mt for this dataset
        mt_path = inputs.as_path(target=sample.dataset, stage=AnnotateDataset, key='mt')
        exp_outputs = self.expected_outputs(sample)

        job = validation_mt_to_vcf_job(
            b=get_batch(),
            mt_path=mt_path,
            sample_id=sample.id,
            out_vcf_path=exp_outputs['vcf'],
            job_attrs=self.get_job_attrs(sample),
            depends_on=inputs.get_jobs(sample),
        )

        return self.make_outputs(sample, data=exp_outputs, jobs=job)
@stage(
    required_stages=ValidationMtToVcf,
    analysis_type='custom',
    update_analysis_meta=_sg_vcf_meta,
    analysis_keys=['vcf'],
)
class ValidationHappyOnVcf(SampleStage):

    def expected_outputs(self, sample: Sample):

        return {
            'vcf': (
                sample.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sample.id}.vcf.bgz'
            ),
            'index': (
                sample.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sample.id}.vcf.bgz.tbi'
            )
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        # only run this on validation samples
        if sample.dataset.name != 'validation':
            return None

        # get the input mt for this dataset
        mt_path = inputs.as_path(target=sample.dataset, stage=AnnotateDataset, key='mt')
        exp_outputs = self.expected_outputs(sample)

        job = validation_mt_to_vcf_job(
            b=get_batch(),
            mt_path=mt_path,
            sample_id=sample.id,
            out_vcf_path=exp_outputs['vcf'],
            job_attrs=self.get_job_attrs(sample),
            depends_on=inputs.get_jobs(sample),
        )

        return self.make_outputs(sample, data=exp_outputs, jobs=job)

