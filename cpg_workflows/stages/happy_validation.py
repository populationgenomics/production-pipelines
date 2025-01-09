#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""

from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.validation import run_happy_on_vcf
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
        job.command(
            f'ss_vcf_from_mt '
            f'--input {input_mt} '
            f'--sample_id {sequencing_group.id} '
            f'--output {str(exp_output)} '
            '--clean ',
        )

        return self.make_outputs(sequencing_group, data=exp_output, jobs=job)


def update_happy_meta(output_path: str) -> dict:
    """
    update
    """
    from csv import DictReader

    from cpg_utils import to_path
    from cpg_utils.config import config_retrieve
    from cpg_workflows.jobs.validation import get_sample_truth_data

    summary_keys = {
        'TRUTH.TOTAL': 'true_variants',
        'METRIC.Recall': 'recall',
        'METRIC.Precision': 'precision',
    }

    happy_handle = to_path(output_path)
    sample_ext_id: str = happy_handle.name.split('__')[0]

    ref_data = get_sample_truth_data(sample_ext_id=sample_ext_id)

    # populate a dictionary of results for this sequencing group
    summary_data = {
        'type': 'validation_result',
        'truth_vcf': ref_data['vcf'],
        'truth_bed': ref_data['bed'],
    }

    if stratification := config_retrieve(['references', 'stratification']):
        summary_data['stratified'] = stratification

    # read in the summary CSV file
    with happy_handle.open() as handle:
        summary_reader = DictReader(handle)
        for line in summary_reader:
            if line['Filter'] != 'PASS' or line['Subtype'] != '*':
                continue

            summary_key = f'{line["Type"]}_{line["Subset"]}'
            for sub_key, sub_value in summary_keys.items():
                summary_data[f'{summary_key}::{sub_value}'] = str(line[sub_key])
    return summary_data


@stage(
    required_stages=ValidationMtToVcf,
    analysis_type='validation',
    analysis_keys=['happy_csv'],
    update_analysis_meta=update_happy_meta,
)
class ValidationHappyOnVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):
        output_prefix = (
            sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / f'{sequencing_group.external_id}__{sequencing_group.id}'
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
        input_vcf = inputs.as_path(target=sequencing_group, stage=ValidationMtToVcf)

        # set the prefix to write outputs to
        output_prefix = (
            sequencing_group.dataset.prefix()
            / 'validation'
            / get_workflow().output_version
            / f'{sequencing_group.external_id}__{sequencing_group.id}'
        )

        exp_outputs = self.expected_outputs(sequencing_group)
        job = run_happy_on_vcf(
            vcf_path=str(input_vcf),
            sample_ext_id=sequencing_group.external_id,
            out_prefix=str(output_prefix),
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=job)
