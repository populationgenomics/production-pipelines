#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.validation import parse_and_post_results, run_happy_on_vcf
from cpg_workflows.stages.talos import query_for_latest_mt
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, get_workflow, stage


@stage()
class ValidationMtToVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        return (
            sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / f'{sequencing_group.id}.vcf.bgz'
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # only keep the sequencing groups with reference data
        if not config_retrieve(['references', sequencing_group.external_id], False):
            return None

        # borrow the talos method to get the latest MT based on analysis entries
        input_mt = config_retrieve(['workflow', 'matrix_table'], query_for_latest_mt(sequencing_group.dataset.name))
        exp_output = self.expected_outputs(sequencing_group)

        job = get_batch().new_job(f'{sequencing_group.id} VCF from dataset MT')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'ss_vcf_from_mt {input_mt} {sequencing_group.id} {str(exp_output)}')

        return self.make_outputs(sequencing_group, data=exp_output, jobs=job)


@stage(required_stages=ValidationMtToVcf, analysis_type='qc', analysis_keys=['happy_csv'])
class ValidationHappyOnVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        output_prefix = (
            sequencing_group.dataset.prefix() / 'validation' / get_workflow().output_version / sequencing_group.id
        )
        return {
            'vcf': to_path(f'{output_prefix}.happy.vcf.bgz'),
            'index': to_path(f'{output_prefix}.happy.vcf.bgz.tbi'),
            'happy_csv': to_path(f'{output_prefix}.happy_extended.csv'),
            'happy_roc': to_path(f'{output_prefix}.happy_roc.all.csv.gz'),
            'happy_metrics': to_path(f'{output_prefix}.happy_metrics.json.gz'),
            'happy_runinfo': to_path(f'{output_prefix}.happy_runinfo.json'),
            'happy_summary': to_path(f'{output_prefix}.summary.csv'),
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # only keep the sequencing groups with reference data
        if not config_retrieve(['references', sequencing_group.external_id], False):
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(target=sequencing_group, stage=ValidationMtToVcf, key='vcf')

        # set the prefix to write outputs to
        output_prefix = (
            sequencing_group.dataset.prefix() / 'validation' / get_workflow().output_version / sequencing_group.id
        )

        exp_outputs = self.expected_outputs(sequencing_group)
        job = run_happy_on_vcf(
            vcf_path=str(input_vcf),
            sequencing_group_ext_id=sequencing_group.external_id,
            out_prefix=str(output_prefix),
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=job)


@stage(
    required_stages=[ValidationMtToVcf, ValidationHappyOnVcf],
    analysis_keys=['json_summary'],
    analysis_type='validation',
)
class ValidationParseHappy(SequencingGroupStage):
    """
    this stage shouldn't exist - we can just parse the results
    of the previous stage using a metadata wrapper
    """

    def expected_outputs(self, sequencing_group: SequencingGroup):
        return {
            'json_summary': sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / f'{sequencing_group.id}.happy_summary.json',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # only keep the sequencing groups with reference data
        if not config_retrieve(['references', sequencing_group.external_id], False):
            return None

        # get the input vcf for this sequence group
        input_vcf = inputs.as_path(target=sequencing_group, stage=ValidationMtToVcf, key='vcf')
        happy_csv = str(inputs.as_dict_by_target(stage=ValidationHappyOnVcf)[sequencing_group.id]['happy_csv'])

        exp_outputs = self.expected_outputs(sequencing_group)

        py_job = get_batch().new_python_job(
            f'parse_{sequencing_group.id}_happy_result',
            (self.get_job_attrs(sequencing_group) or {}) | {'tool': 'hap.py'},
        )
        py_job.image(get_config()['workflow']['driver_image'])
        py_job.call(
            parse_and_post_results,
            str(input_vcf),
            sequencing_group.id,
            sequencing_group.external_id,
            happy_csv,
            str(exp_outputs['json_summary']),
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=py_job)
