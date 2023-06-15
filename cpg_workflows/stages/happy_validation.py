#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""

import logging
from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    SequencingGroupStage,
    SequencingGroup,
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
class ValidationMtToVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {
            'vcf': (
                sequencing_group.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sequencing_group.id}.vcf.bgz'
            ),
            'index': (
                sequencing_group.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / f'{sequencing_group.id}.vcf.bgz.tbi'
            ),
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # only run this on validation sequencing groups
        if sequencing_group.dataset.name != 'validation':
            return None

        # only keep the sequencing groups with reference data
        if sequencing_group.external_id not in get_config()['references']:
            return None

        # get the input mt for this dataset
        mt_path = inputs.as_path(
            target=sequencing_group.dataset, stage=AnnotateDataset, key='mt'
        )
        exp_outputs = self.expected_outputs(sequencing_group)

        job = validation_mt_to_vcf_job(
            b=get_batch(),
            mt_path=str(mt_path),
            sequencing_group_id=sequencing_group.id,
            out_vcf_path=str(exp_outputs['vcf']),
            job_attrs=self.get_job_attrs(sequencing_group),
            depends_on=inputs.get_jobs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=job)


@stage(
    required_stages=ValidationMtToVcf,
    analysis_type='qc',
    analysis_keys=['happy_csv'],
)
class ValidationHappyOnVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        output_prefix = (
            sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / sequencing_group.id
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

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # only run this on validation sequencing groups
        if sequencing_group.dataset.name != 'validation':
            return None

        # only keep the sequencing groups with reference data
        if sequencing_group.external_id not in get_config()['references']:
            logging.info(f'Skipping {sequencing_group.id}; not in the reference set')
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(
            target=sequencing_group, stage=ValidationMtToVcf, key='vcf'
        )

        # set the prefix to write outputs to
        output_prefix = (
            sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / sequencing_group.id
        )

        exp_outputs = self.expected_outputs(sequencing_group)
        job = run_happy_on_vcf(
            b=get_batch(),
            vcf_path=str(input_vcf),
            sequencing_group_ext_id=sequencing_group.external_id,
            out_prefix=str(output_prefix),
            job_attrs=self.get_job_attrs(sequencing_group),
            depends_on=inputs.get_jobs(sequencing_group),
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=job)


@stage(required_stages=[ValidationMtToVcf, ValidationHappyOnVcf])
class ValidationParseHappy(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {
            'json_summary': (
                sequencing_group.dataset.prefix()
                / 'validation'
                / get_workflow().output_version
                / sequencing_group.id,
                'happy_summary.json',
            )
        }

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        # only run this on validation sequencing groups
        if sequencing_group.dataset.name != 'validation':
            return None

        # only keep the sequencing groups with reference data
        if sequencing_group.external_id not in get_config()['references']:
            logging.info(f'Skipping {sequencing_group.id}; not in the reference set')
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(
            target=sequencing_group, stage=ValidationMtToVcf, key='vcf'
        )
        happy_results = inputs.as_dict_by_target(stage=ValidationHappyOnVcf)[
            sequencing_group.id
        ]

        exp_outputs = self.expected_outputs(sequencing_group)

        py_job = get_batch().new_python_job(
          f'parse_{sequencing_group.id}_happy_result',
          (self.get_job_attrs(sequencing_group) or {}) | {'tool': 'hap.py'}
        )
        py_job.image(get_config()['workflow']['driver_image'])
        py_job.call(
            parse_and_post_results,
            vcf_path=str(input_vcf),
            sequencing_group_id=sequencing_group.external_id,
            happy_results=str(happy_results),
            out_file=str(exp_outputs['json_summary'])
        )
        # set dependencies if applicable
        if dependencies := inputs.get_jobs(sequencing_group):
            py_job.depends_on(*dependencies)

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=py_job)
