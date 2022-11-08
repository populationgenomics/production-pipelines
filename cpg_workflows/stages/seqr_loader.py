#!/usr/bin/env python3

"""
Hail Query stages for the Seqr loader workflow.
"""

from cpg_utils import to_path, Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
    Cohort,
    Dataset,
    get_workflow,
)
from cpg_workflows.jobs.seqr_loader import annotate_cohort_jobs, annotate_dataset_jobs

from .joint_genotyping import JointGenotyping
from .vep import Vep
from .vqsr import Vqsr
from .. import get_cohort, get_batch


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
            'tmp_prefix': str(self.tmp_prefix / 'mt' / h),
            'mt': get_cohort().analysis_dataset.prefix() / 'mt' / f'{h}.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, key='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(
            target=cohort, stage=Vqsr, key='siteonly'
        )
        vep_ht_path = inputs.as_path(target=cohort, stage=Vep, key='ht')

        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['tmp_prefix']) / 'checkpoints'
        )

        jobs = annotate_cohort_jobs(
            b=get_batch(),
            vcf_path=vcf_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_vcf_path,
            vep_ht_path=vep_ht_path,
            out_mt_path=self.expected_outputs(cohort)['mt'],
            checkpoint_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(cohort),
            depends_on=inputs.get_jobs(cohort),
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
        h = get_cohort().alignment_inputs_hash()
        return {
            'tmp_prefix': str(self.tmp_prefix / 'mt' / f'{h}-{dataset.name}'),
            'mt': dataset.prefix() / 'mt' / f'{h}-{dataset.name}.mt',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        assert dataset.cohort
        mt_path = inputs.as_path(target=dataset.cohort, stage=AnnotateCohort, key='mt')

        checkpoint_prefix = (
            to_path(self.expected_outputs(dataset)['tmp_prefix']) / 'checkpoints'
        )

        jobs = annotate_dataset_jobs(
            b=get_batch(),
            mt_path=mt_path,
            sample_ids=dataset.get_sample_ids(),
            out_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset),
        )

        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


def es_password() -> str:
    """
    Get Elasticsearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return read_secret(
        project_id=get_config()['elasticsearch']['password_project_id'],
        secret_name=get_config()['elasticsearch']['password_secret_id'],
        fail_gracefully=False,
    )


@stage(
    required_stages=[AnnotateDataset],
    analysis_type='es-index',
    analysis_key='index_name',
)
class MtToEs(DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = get_config()['workflow']['sequencing_type']
        index_name = (
            f'{dataset.name}-{sequencing_type}-{get_workflow().run_timestamp}'.lower()
        )
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        if (
            es_datasets := get_config()['workflow'].get('create_es_index_for_datasets')
        ) and dataset.name not in es_datasets:
            # Skipping dataset that wasn't explicitly requested to upload to ES:
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(
            target=dataset, stage=AnnotateDataset, key='mt'
        )
        index_name = self.expected_outputs(dataset)['index_name']
        done_flag_path = self.expected_outputs(dataset)['done_flag']

        if 'elasticsearch' not in get_config():
            raise ValueError(
                f'"elasticsearch" section is not defined in config, cannot create '
                f'Elasticsearch index for dataset {dataset}'
            )

        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        script_path = to_path(__file__).parent / 'dataproc_scripts' / 'mt_to_es.py'

        j = dataproc.hail_dataproc_job(
            get_batch(),
            f'{script_path} '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {index_name} '
            f'--done-flag-path {done_flag_path} '
            f'--es-password {es_password()} '
            f'--liftover-path {reference_path("liftover_38_to_37")}',
            max_age='24h',
            packages=[
                'cpg_workflows',
                'elasticsearch==8.*',
                'google',
                'fsspec',
                'gcloud',
            ],
            num_workers=2,
            num_secondary_workers=0,
            job_name=f'{dataset.name}: create ES index',
            depends_on=inputs.get_jobs(dataset),
            scopes=['cloud-platform'],
            pyfiles=['seqr-loading-pipelines/hail_scripts'],
        )
        j._preemptible = False
        j.attributes = (j.attributes or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]
        return self.make_outputs(dataset, data=index_name, jobs=jobs)
