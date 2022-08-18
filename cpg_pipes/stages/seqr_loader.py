#!/usr/bin/env python3

"""
Seqr loading pipeline: FASTQ -> ElasticSearch index.
"""

import logging

from cpg_utils import to_path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config

from cpg_pipes import utils
from cpg_pipes.jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from cpg_pipes.pipeline import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
)
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vep import Vep
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.targets import Cohort, Dataset
from cpg_utils.hail_batch import reference_path

logger = logging.getLogger(__file__)


@stage(required_stages=[JointGenotyping, Vqsr, Vep])
class AnnotateCohort(CohortStage):
    """
    Re-annotate the entire cohort.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a matrix table.
        """
        h = cohort.alignment_inputs_hash()
        return {
            'prefix': str(self.tmp_prefix / 'mt' / h),
            'mt': cohort.analysis_dataset.tmp_prefix() / 'mt' / f'{h}.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, id='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(
            target=cohort, stage=Vqsr, id='siteonly'
        )
        vep_ht_path = inputs.as_path(target=cohort, stage=Vep, id='ht')

        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['prefix']) / 'checkpoints'
        )

        jobs = annotate_cohort_jobs(
            b=self.b,
            vcf_path=vcf_path,
            vep_ht_path=vep_ht_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_vcf_path,
            output_mt_path=self.expected_outputs(cohort)['mt'],
            checkpoint_prefix=checkpoint_prefix,
            sequencing_type=get_config()['workflow']['sequencing_type'],
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(
            cohort,
            data=self.expected_outputs(cohort),
            jobs=jobs,
        )


@stage(required_stages=[AnnotateCohort])
class AnnotateDataset(DatasetStage):
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """

    def expected_outputs(self, dataset: Dataset):
        """
        Expected to generate a matrix table
        """
        h = self.cohort.alignment_inputs_hash()
        return {
            'prefix': str(self.tmp_prefix / 'mt' / f'{h}-{dataset.name}'),
            # We want to write the matrix table into the main bucket.
            'mt': dataset.prefix() / 'mt' / f'{h}-{dataset.name}.mt',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        mt_path = inputs.as_path(target=self.cohort, stage=AnnotateCohort, id='mt')

        checkpoint_prefix = (
            to_path(self.expected_outputs(dataset)['prefix']) / 'checkpoints'
        )

        jobs = annotate_dataset_jobs(
            b=self.b,
            mt_path=mt_path,
            sample_ids=[s.id for s in dataset.get_samples()],
            output_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_bucket=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
            overwrite=not get_config()['workflow'].get('check_intermediates'),
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


def es_password() -> str:
    """
    Get ElasticSearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return read_secret(
        project_id=get_config()['elasticsearch']['password_project_id'],
        secret_name=get_config()['elasticsearch']['password_secret_id'],
        fail_gracefully=False,
    )


@stage(required_stages=[AnnotateDataset], analysis_type='es-index')
class LoadToEs(DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: Dataset) -> str:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = get_config()['workflow']['sequencing_type']
        index_name = f'{dataset.name}-{sequencing_type}-{self.run_id}'
        return index_name

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        if (
            es_datasets := get_config()['workflow'].get('create_es_index_for_datasets')
        ) and dataset.name not in es_datasets:
            # Skipping dataset that wasn't explicitly requested to upload to ES:
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDataset, id='mt')
        index_name = self.expected_outputs(dataset).lower()

        from analysis_runner import dataproc

        j = dataproc.hail_dataproc_job(
            self.b,
            f'cpg_pipes/dataproc_scripts/seqr/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {index_name} '
            f'--es-password {es_password()} '
            f'--liftover-path {reference_path("liftover_38_to_37")} '
            f'--use-spark ',  # es export doesn't work with the service backend
            max_age='24h',
            packages=utils.DATAPROC_PACKAGES,
            num_workers=2,
            num_secondary_workers=0,
            job_name=f'{dataset.name}: create ES index',
            depends_on=inputs.get_jobs(dataset),
            scopes=['cloud-platform'],
        )
        j._preemptible = False
        j.attributes = self.get_job_attrs(dataset) | {'tool': 'hail dataproc'}
        jobs = [j]
        return self.make_outputs(dataset, data=index_name, jobs=jobs)
