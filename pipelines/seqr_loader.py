#!/usr/bin/env python3

"""
Seqr loading pipeline: FASTQ -> ElasticSearch index.
"""

import logging
import os
import time

import click
import yaml
from google.cloud import secretmanager

from cpg_utils.config import get_config

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from cpg_pipes.pipeline import (
    pipeline_entry_point,
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
)
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.vep import Vep
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.targets import Cohort, Dataset

logger = logging.getLogger(__file__)

SUPPORTED_SEQUENCING_TYPES = ['genome']


@stage(required_stages=[JointGenotyping, Vep, Vqsr])
class AnnotateCohort(CohortStage):
    """
    Re-annotate the entire cohort.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expected to write a matrix table.
        """
        h = cohort.alignment_inputs_hash()
        return cohort.analysis_dataset.path() / 'mt' / f'{h}.mt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, id='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(target=cohort, stage=Vqsr)
        vep_ht_path = inputs.as_path(target=cohort, stage=Vep)

        mt_path = self.expected_outputs(cohort)

        jobs = annotate_cohort_jobs(
            b=self.b,
            vcf_path=vcf_path,
            vep_ht_path=vep_ht_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_vcf_path,
            output_mt_path=mt_path,
            checkpoints_bucket=self.tmp_bucket / 'checkpoints',
            sequencing_type=cohort.get_sequencing_type(),
            overwrite=not get_config()['workflow'].get('self.check_intermediates'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=mt_path, jobs=jobs)


@stage(required_stages=[AnnotateCohort])
class AnnotateDataset(DatasetStage):
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        h = self.cohort.alignment_inputs_hash()
        return self.tmp_bucket / f'{h}-{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        mt_path = inputs.as_path(target=self.cohort, stage=AnnotateCohort)

        jobs = annotate_dataset_jobs(
            b=self.b,
            mt_path=mt_path,
            sample_ids=[s.id for s in dataset.get_samples()],
            output_mt_path=self.expected_outputs(dataset),
            tmp_bucket=self.tmp_bucket / 'checkpoints' / dataset.name,
            job_attrs=self.get_job_attrs(dataset),
            overwrite=not get_config()['workflow'].get('self.check_intermediates'),
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(required_stages=[AnnotateDataset])
class LoadToEs(DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: Dataset) -> None:
        """
        Expected to generate a Seqr index, which is not a file
        """
        return None

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        if (
            es_datasets := get_config()['workflow'].get('create_es_index_for_datasets')
        ) and dataset.name not in es_datasets:
            # Skipping dataset that wasn't explicitly requested to upload to ES:
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDataset)
        version = time.strftime('%Y%m%d-%H%M%S')
        index_name = f'{dataset.name}-{version}'
        
        from analysis_runner import dataproc
        j = dataproc.hail_dataproc_job(
            self.b,
            f'cpg_pipes/dataproc_scripts/seqr/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-host elasticsearch.es.australia-southeast1.gcp.elastic-cloud.com '
            f'--es-port 9243 '
            f'--es-username seqr '
            f'--es-password {_read_es_password()} '
            f'--es-index {index_name} '
            f'--es-index-min-num-shards 1 '
            f'--use-spark ',  # es export doesn't work with the service backend
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=2,
            job_name=f'{dataset.name}: create ES index',
            depends_on=inputs.get_jobs(dataset),
            scopes=['cloud-platform'],
        )
        jobs = [j]
        if self.status_reporter:
            jobs = self.status_reporter.add_updaters_jobs(
                self.b,
                output=index_name,
                analysis_type='seqr_index',
                target=dataset,
                jobs=jobs,
            )
        j.attributes = self.get_job_attrs(dataset)
        return self.make_outputs(dataset, jobs=jobs)


def _read_es_password(
    project_id='seqr-308602',
    secret_id='seqr-es-password',
) -> str:
    """
    Read a GCP secret storing the ES password
    """
    if password := os.environ.get('SEQR_ES_PASSWORD'):
        return password
    client = secretmanager.SecretManagerServiceClient()
    secret_path = client.secret_version_path(project_id, secret_id, 'latest')
    # noinspection PyTypeChecker
    response = client.access_secret_version(request={'name': secret_path})
    return response.payload.data.decode('UTF-8')


@pipeline_entry_point(name='Seqr Loader')
def main(pipeline: Pipeline):
    """
    Seqr loading pipeline: FASTQ -> ElasticSearch index.
    """
    if seq_type := get_config()['workflow'].get('sequencing_type'):
        if seq_type not in SUPPORTED_SEQUENCING_TYPES:
            raise click.BadParameter(
                f'Unsupported sequencing data type {seq_type.value}. '
                f'Supported types: {[st for st in SUPPORTED_SEQUENCING_TYPES]} '
            )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
