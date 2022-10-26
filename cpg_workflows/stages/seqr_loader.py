#!/usr/bin/env python3

"""
Hail Query stages for the Seqr loader workflow.
"""

from cpg_utils import to_path, Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path
from cpg_utils.workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
    Cohort,
    Dataset,
    get_workflow,
)

from .joint_genotyping import JointGenotyping
from .vep import Vep
from .vqsr import Vqsr


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
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, key='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(
            target=cohort, stage=Vqsr, key='siteonly'
        )
        vep_ht_path = inputs.as_path(target=cohort, stage=Vep, key='ht')

        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['prefix']) / 'checkpoints'
        )

        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        j = dataproc.hail_dataproc_job(
            self.b,
            f'dataproc_scripts/annotate_cohort.py '
            f'--vcf-path {vcf_path} '
            f'--siteonly-vqsr-vcf-path {siteonly_vqsr_vcf_path} '
            f'--vep-ht-path {vep_ht_path} '
            f'--out-mt-path {self.expected_outputs(cohort)["mt"]} '
            f'--checkpoint-prefix {checkpoint_prefix}',
            max_age='24h',
            packages=[
                'cpg_utils',
                'google',
                'fsspec',
                'gcloud',
            ],
            num_workers=2,
            num_secondary_workers=20,
            job_name=f'Annotate cohort',
            depends_on=inputs.get_jobs(cohort),
            scopes=['cloud-platform'],
            pyfiles=[
                'seqr-loading-pipelines/hail_scripts',
                'query_modules',
            ],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        )
        j.attributes = self.get_job_attrs(cohort) | {'tool': 'hailctl dataproc'}

        return self.make_outputs(
            cohort,
            data=self.expected_outputs(cohort),
            jobs=[j],
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
        mt_path = inputs.as_path(target=self.cohort, stage=AnnotateCohort, key='mt')

        checkpoint_prefix = (
            to_path(self.expected_outputs(dataset)['prefix']) / 'checkpoints'
        )

        sample_ids_list_path = dataset.tmp_prefix() / 'sample-list.txt'
        if not get_config()['hail'].get('dry_run', False):
            with sample_ids_list_path.open('w') as f:
                f.write(','.join([s.id for s in dataset.get_samples()]))

        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        j = dataproc.hail_dataproc_job(
            self.b,
            f'dataproc_scripts/annotate_dataset.py '
            f'--mt-path {mt_path} '
            f'--sample-ids {sample_ids_list_path} '
            f'--out-mt-path {self.expected_outputs(dataset)["mt"]} '
            f'--checkpoint-prefix {checkpoint_prefix}',
            max_age='24h',
            packages=[
                'cpg_utils',
                'google',
                'fsspec',
                'gcloud',
            ],
            num_workers=2,
            num_secondary_workers=20,
            job_name=f'Annotate dataset',
            depends_on=inputs.get_jobs(dataset),
            scopes=['cloud-platform'],
            pyfiles=[
                'seqr-loading-pipelines/hail_scripts',
                'query_modules',
            ],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
        )
        j.attributes = self.get_job_attrs(dataset) | {'tool': 'hailctl dataproc'}

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[j])


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


@stage(required_stages=[AnnotateDataset], analysis_type='es-index')
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

        j = dataproc.hail_dataproc_job(
            self.b,
            f'dataproc_scripts/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {index_name} '
            f'--done-flag-path {done_flag_path} '
            f'--es-password {es_password()} '
            f'--liftover-path {reference_path("liftover_38_to_37")}',
            max_age='24h',
            packages=[
                'cpg_utils',
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
