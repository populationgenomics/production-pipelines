"""
Driver for computing structural variants from gatk-sv.
Uses metadata from the sample-metadata server.

- 2021-10-28 Michael Franklin and Vlad Savelyev
"""
from collections import defaultdict
from typing import Collection, Optional, List, Dict

import os

import click
import hailtop.batch as hb
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)

from .util import SMDB


SV_WORKFLOWS_COMMIT = ''

STAGES = [
    # put stages here of pipeline
]

VARIANT_TYPES = ['manta', 'melt', 'wham']

DATASET = os.getenv('DATASET')
BUCKET = os.getenv('HAIL_BUCKET')
OUTPUT_SUFFIX = os.getenv('OUTPUT')
BILLING_PROJECT = os.getenv('HAIL_BILLING_PROJECT')
ACCESS_LEVEL = os.getenv('ACCESS_LEVEL')


@click.command()
def main(
    analysis_dataset: str,
    input_datasets: Collection[str],
    start_from_stage: Optional[str],
    end_with_stage: Optional[str],
    namespace: str,
    skip_samples: Collection[str],
):
    assert input_datasets

    _skip_samples = set(skip_samples) if skip_samples else set()

    if namespace == 'test':
        sample_metadata_projects = [d + '-test' for d in input_datasets]
    else:
        sample_metadata_projects = input_datasets

    _samples_by_dataset = {}
    _sample_ids_by_dataset = {}
    for dataset, sm_project in zip(input_datasets, sample_metadata_projects):
        samples = SMDB.get_samples_by_project(
            project=sm_project, skip_samples=_skip_samples
        )
        _samples_by_dataset[dataset] = samples
        _samples_by_dataset[dataset] = [s['id'] for s in samples]

    b = hb.Batch('gatk-sv-' + '-'.join(input_datasets))
    # gather batch evidence
    run_gather_batch_evidence(
        b,
        analysis_dataset=analysis_dataset,
        input_datasets=input_datasets,
        sample_metadata_projects=sample_metadata_projects,
        sample_ids_by_dataset=_sample_ids_by_dataset,
    )


def run_gather_batch_evidence(
    b: hb.Batch,
    analysis_dataset: str,
    input_datasets: Collection[str],
    sample_metadata_projects: Collection[str],
    sample_ids_by_dataset: Dict[str, List[str]],
):
    _crams_by_sample_id = {}
    _variants_by_sample_type = defaultdict(dict)
    for dataset, sm_project in zip(input_datasets, sample_metadata_projects):
        sample_ids = sample_ids_by_dataset[dataset]
        for variant_type in VARIANT_TYPES:
            analyses = SMDB.get_variants_for_sample_id_by_sv_type(
                project=sm_project, sample_ids=sample_ids, sv_type=variant_type
            )
            for a in analyses:
                sample_id = a['sample_ids'][0]
                _variants_by_sample_type[sample_id][variant_type] = a

        cram_analyses = SMDB.get_cram_analyses_from_project_for_sample_ids(
            project=sm_project, sample_ids=sample_ids
        )
        for a in cram_analyses:
            sample_id = a['sample_ids'][0]
            _crams_by_sample_id[sample_id] = a

    for dataset, sm_project in zip(input_datasets, sample_metadata_projects):
        _samples_ids = sample_ids_by_dataset[dataset]
        # this is tricky, so we'll actually just grab ALL samples that don't have
        # ALL of the associated variant_types
        _sample_ids_for_gather_batch_evidence = [
            sid
            for sid in _samples_ids
            if not all(
                variant_type in _variants_by_sample_type[sid]
                for variant_type in VARIANT_TYPES
            )
        ]

        cram_paths = [
            _crams_by_sample_id[sid]['output']
            for sid in _sample_ids_for_gather_batch_evidence
        ]

        nouts = len(_sample_ids_for_gather_batch_evidence)
        gather_batched_evidence = run_cromwell_workflow_from_repo_and_get_outputs(
            b=b,
            job_prefix='gather_batch_evidence',
            dataset=analysis_dataset,
            access_level=ACCESS_LEVEL,
            repo='gatk-sv',
            commit=SV_WORKFLOWS_COMMIT,
            cwd='gatk-sv/test/gatk-sv-git/wdl',
            workflow='GatherSampleEvidenceBatch.wdl',
            libs=['.'],
            output_suffix=OUTPUT_SUFFIX,
            # this has to be relative to the sv-workflows repository
            input_paths=[],
            input_dict={
                'GatherSampleEvidenceBatch.sample_ids': _sample_ids_for_gather_batch_evidence,
                'GatherSampleEvidenceBatch.bam_or_cram_files': cram_paths,
            },
            outputs_to_collect={
                # output paths
                'manta_vcf_paths': CromwellOutputType.array_path(
                    name='GatherSampleEvidenceBatch.manta_vcf', length=nouts
                ),
                'melt_vcf_paths': CromwellOutputType.array_path(
                    name='GatherSampleEvidenceBatch.melt_vcf', length=nouts
                ),
                'wham_vcf_paths': CromwellOutputType.array_path(
                    name='GatherSampleEvidenceBatch.wham_vcf', length=nouts
                ),
                # output files
                # 'coverage_counts': CromwellOutputType.array(
                #     name='GatherSampleEvidenceBatch.coverage_counts', length=nouts
                # ),
                # 'manta_vcf': CromwellOutputType.array_resource_group(
                #     name='GatherSampleEvidenceBatch.manta_vcf',
                #     length=nouts,
                #     resource_group={
                #         'vcf.gz': 'GatherSampleEvidenceBatch.manta_vcf',
                #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.manta_index',
                #     },
                # ),
                # 'melt_vcf': CromwellOutputType.array_resource_group(
                #     name='GatherSampleEvidenceBatch.melt_vcf',
                #     length=nouts,
                #     resource_group={
                #         'vcf.gz': 'GatherSampleEvidenceBatch.melt_vcf',
                #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.melt_index',
                #     },
                # ),
                # 'wham_vcf': CromwellOutputType.array_resource_group(
                #     name='GatherSampleEvidenceBatch.wham_vcf',
                #     length=nouts,
                #     resource_group={
                #         'vcf.gz': 'GatherSampleEvidenceBatch.wham_vcf',
                #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.wham_index',
                #     },
                # ),
                # 'melt_coverage': CromwellOutputType.array(
                #     name='GatherSampleEvidenceBatch.melt_coverage', length=nouts
                # ),
                # 'melt_read_length': CromwellOutputType.array(
                #     name='GatherSampleEvidenceBatch.melt_read_length', length=nouts
                # ),
                # 'melt_insert_size': CromwellOutputType.array(
                #     name='GatherSampleEvidenceBatch.melt_insert_size', length=nouts
                # ),
                # 'pesr_disc': CromwellOutputType.array_resource_group(
                #     name='GatherSampleEvidenceBatch.pesr_disc',
                #     length=nouts,
                #     resource_group={
                #         'txt.gz': 'GatherSampleEvidenceBatch.pesr_disc',
                #         'txt.gz.tbi': 'GatherSampleEvidenceBatch.pesr_disc_index',
                #     },
                # ),
                # 'pesr_split': CromwellOutputType.array_resource_group(
                #     name='GatherSampleEvidenceBatch.pesr_split',
                #     length=nouts,
                #     resource_group={
                #         'txt.gz': 'GatherSampleEvidenceBatch.pesr_split',
                #         'txt.gz.tbi': 'GatherSampleEvidenceBatch.pesr_split_index',
                #     },
                # ),
                # # optional
                # 'metrics_file_sampleevidence': CromwellOutputType.single(
                #     name='GatherSampleEvidenceBatch.metrics_file_sampleevidence',
                # ),
            },
        )

        create_sm_sv_analyis_object_from_batch(
            b,
            'manta',
            sm_project,
            _sample_ids_for_gather_batch_evidence,
            gather_batched_evidence['manta_vcf_paths'],
        )
        create_sm_sv_analyis_object_from_batch(
            b,
            'melt',
            sm_project,
            _sample_ids_for_gather_batch_evidence,
            gather_batched_evidence['melt_vcf_paths'],
        )
        create_sm_sv_analyis_object_from_batch(
            b,
            'wham',
            sm_project,
            _sample_ids_for_gather_batch_evidence,
            gather_batched_evidence['wham_vcf_paths'],
        )


# this is a little tricky, because this will actually pass the file to it ://
def create_sm_sv_analyis_object_from_batch(
    b: hb.Batch,
    sm_project: str,
    sv_type: str,
    sample_ids_: List[str],
    file_resources_containing_path: List,
):
    """
    Create a job (within a batch) that creates analysis objects (of type=sv)
    for each pair of (sample_id, file_resource_containing_path), where the
    file_resource_containing_path might be returned from Cromwell
    """

    def sm_create_sv_analysis_call(sample_id_: str, file_containing_path):
        """
        Batch callable function that creates SV analysis object on SM server
        """
        with open(file_containing_path) as f:
            vcf_path = file_containing_path.read().strip()

        # consider moving the file before creating_analysis, eg:
        #   permanent_location = 'gs://<bucket>/{sample_id}.sv.{sv_type}.vcf.gz'
        #   gsutil mv {vcf_path} {permanent_location}
        SMDB.create_analysis(
            sm_project,
            'sv',
            vcf_path,
            'completed',
            [sample_id_],
            {'sv_algorithm': sv_type},
        )

    j = b.new_python_job(f'sm-update-{sv_type}-path')
    for sample_id, file_resource_containing_path in zip(
        sample_ids_, file_resources_containing_path
    ):
        j.call(
            sm_create_sv_analysis_call,
            sm_project,
            sample_id,
            sv_type,
            file_resource_containing_path,
        )
