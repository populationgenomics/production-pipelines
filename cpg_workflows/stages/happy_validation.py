#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""


from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    SampleStage,
    Sample,
    StageInput,
    StageOutput,
    get_workflow,
)
from cpg_workflows.stages.seqr_loader import AnnotateDataset, _sg_vcf_meta
from cpg_workflows.jobs.validation import (
    validation_mt_to_vcf_job,
    run_happy_on_vcf,
    parse_and_post_results,
)
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
            ),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        # only run this on validation samples
        if sample.dataset.name != 'validation':
            return None

        # only keep the samples with reference data
        if sample.external_id not in get_config()['references']:
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
    analysis_type='qc',
    analysis_keys=['happy_csv'],
)
class ValidationHappyOnVcf(SampleStage):
    def expected_outputs(self, sample: Sample):
        output_prefix = (
            sample.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / sample.id
        )
        return {
            'vcf': output_prefix / 'happy.vcf.bgz',
            'index': output_prefix / 'happy.vcf.bgz.tbi',
            'happy_csv': output_prefix / 'happy_extended.csv',
            'happy_roc': output_prefix / 'happy_roc.all.csv.gz',
            'happy_metrics': output_prefix / 'happy_metrics.json.gz',
            'happy_runinfo': output_prefix / 'happy_runinfo.json',
            'happy_summary': output_prefix / 'summary.csv',
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        # only run this on validation samples
        if sample.dataset.name != 'validation':
            return None

        # only keep the samples with reference data
        if sample.external_id not in get_config()['references']:
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(target=sample, stage=ValidationMtToVcf, key='vcf')

        # set the prefix to write outputs to
        output_prefix = (
            sample.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / sample.id
        )

        exp_outputs = self.expected_outputs(sample)
        job = run_happy_on_vcf(
            b=get_batch(),
            vcf_path=str(input_vcf),
            sample_ext_id=sample.external_id,
            out_prefix=str(output_prefix),
            job_attrs=self.get_job_attrs(sample),
            depends_on=inputs.get_jobs(sample),
        )

        return self.make_outputs(sample, data=exp_outputs, jobs=job)


@stage(required_stages=[ValidationMtToVcf, ValidationHappyOnVcf])
class ValidationParseHappy(SampleStage):
    def expected_outputs(self, sample: Sample):
        return {
            'json_summary': (
                sample.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / sample.id,
                'happy_summary.json',
            )
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        # only run this on validation samples
        if sample.dataset.name != 'validation':
            return None

        # only keep the samples with reference data
        if sample.external_id not in get_config()['references']:
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(target=sample, stage=ValidationMtToVcf, key='vcf')
        happy_results = inputs.as_dict_by_target(stage=ValidationHappyOnVcf)[sample.id]

        exp_outputs = self.expected_outputs(sample)
        job = parse_and_post_results(
            b=get_batch(),
            vcf_path=str(input_vcf),
            sample=sample,
            happy_results=happy_results,
            out_file=exp_outputs['json_summary'],
            job_attrs=self.get_job_attrs(sample),
            depends_on=inputs.get_jobs(sample),
        )

        return self.make_outputs(sample, data=exp_outputs, jobs=job)
