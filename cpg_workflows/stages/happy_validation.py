#!/usr/bin/env python3

"""
Content relating to the hap.py validation process
"""
import re
from functools import lru_cache

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.validation import run_happy_on_vcf
from cpg_workflows.stages.talos import query_for_latest_hail_object
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage

# this pipeline is only ever expected to run for one project
VALIDATION_PROJECT: str = 'validation'
# default Stage-as-source is the RD_Combiner version
ANNOTATE_DATASET: str = 'AnnotateDatasetSmallVariantsWithHailQuery'
# a regex to pull the output_version hash from the full output path
OUTPUT_VERSION_REGEX = re.compile(r'gs://.+/(?P<hash>\w+_\d+)')


def find_hash_from_path(input_path: str) -> str | ValueError:
    """
    Extracts the Hash from a Google cloud Path

    Args:
        input_path (str): full path to a file/folder in GCP

    Returns:
        A string if it was successfully located, or a ValueError Exception type
    """

    result = re.match(OUTPUT_VERSION_REGEX, input_path)

    if not result:
        raise ValueError(f'Unable to identify a Hash in {input_path}')

    return result['hash']


@lru_cache(1)
def find_input_mt_and_hash():
    """
    Pulls the input MT to use for this project from Metamist (or config)
    Returns both the path and the hash component of it

    Returns:
        2 strings if successful, or a ValueError
    """

    # if there's a specific override in config, use that
    if input_mt := config_retrieve(['workflow', 'matrixtable'], None):
        # get a chunk of the hash from the String
        return input_mt, find_hash_from_path(input_mt)

    # borrow the talos method to get the latest MT based on analysis entries
    # we allow overriding from config to use a specific table
    # otherwise we query for the latest MT, using analysis_type to allow swapping between rd_combiner & seqr_loader
    input_mt = query_for_latest_hail_object(
        VALIDATION_PROJECT,
        analysis_type=config_retrieve(['workflow', 'mt_entry_type'], default='matrixtable'),
        object_suffix='.mt',
        exact_string=config_retrieve(['workflow', 'mt_stage'], ANNOTATE_DATASET),
    )

    return input_mt, find_hash_from_path(input_mt)


@stage()
class ValidationMtToVcf(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup):

        _mt, hash_str = find_input_mt_and_hash()

        return to_path(
            dataset_path(
                dataset=VALIDATION_PROJECT,
                # repeat "validation" as pipeline folder, within the validation bucket
                suffix=f'validation/{hash_str}/{sequencing_group.id}.vcf.bgz',
            ),
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        # only keep the sequencing groups with reference data
        if not config_retrieve(['references', sequencing_group.external_id], False):
            return None

        input_mt, hash_str = find_input_mt_and_hash()

        exp_output = self.expected_outputs(sequencing_group)

        job = get_batch().new_job(f'{sequencing_group.id} VCF from MT ({hash_str})')
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
        _mt, hash_str = find_input_mt_and_hash()

        output_prefix = dataset_path(
            dataset=VALIDATION_PROJECT,
            # repeat "validation" as pipeline folder, within the validation bucket
            suffix=f'validation/{hash_str}/{sequencing_group.external_id}__{sequencing_group.id}',
        )

        return {
            'prefix': output_prefix,
            'vcf': to_path(f'{output_prefix}.happy.vcf.gz'),
            'index': to_path(f'{output_prefix}.happy.vcf.gz.tbi'),
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

        exp_outputs = self.expected_outputs(sequencing_group)

        job = run_happy_on_vcf(
            vcf_path=str(input_vcf),
            sample_ext_id=sequencing_group.external_id,
            out_prefix=exp_outputs['prefix'],
        )

        return self.make_outputs(sequencing_group, data=exp_outputs, jobs=job)
