#!/usr/bin/env python3

"""
Hail Query stages for the Seqr loader workflow.
"""
from typing import Any

from cpg_utils import to_path, Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config
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
from cpg_workflows.jobs.seqr_loader import (
    annotate_cohort_jobs,
    annotate_dataset_jobs,
    cohort_to_vcf_job,
)

from .joint_genotyping import JointGenotyping
from .vep import Vep
from .vqsr import Vqsr
from .. import get_batch


@stage(required_stages=[JointGenotyping, Vqsr, Vep])
class AnnotateCohort(CohortStage):
    """
    Re-annotate the entire cohort.
    """

    def expected_outputs(self, cohort: Cohort):
        """
        Expected to write a matrix table.
        """
        return {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.tmp_prefix),
            'mt': self.prefix / 'cohort.mt',
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


def _update_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to custom analysis object

    TODO: Replace this once dynamic analysis types land in metamist.
    """
    return {'type': 'annotated-dataset-callset'}


def _dataset_vcf_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to custom analysis object

    TODO: Replace this once dynamic analysis types land in metamist.
    """
    return {'type': 'dataset-vcf'}


def _sg_vcf_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to custom analysis object

    TODO: Replace this once dynamic analysis types land in metamist.
    """
    return {'type': 'dataset-vcf'}


@stage(
    required_stages=[AnnotateCohort],
    analysis_type='custom',
    update_analysis_meta=_update_meta,
    analysis_keys=['mt'],
)
class AnnotateDataset(DatasetStage):
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """

    def expected_outputs(self, dataset: Dataset):
        """
        Expected to generate a matrix table
        """
        return {
            'tmp_prefix': str(self.tmp_prefix / dataset.name),
            'mt': (
                dataset.prefix()
                / 'mt'
                / f'{get_workflow().output_version}-{dataset.name}.mt'
            ),
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
            sequencing_group_ids=dataset.get_sequencing_group_ids(),
            out_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset),
        )

        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(
    required_stages=[AnnotateDataset],
    analysis_type='custom',
    update_analysis_meta=_dataset_vcf_meta,
    analysis_keys=['vcf'],
)
class DatasetVCF(DatasetStage):
    """
    Take the per-dataset MT and write out as a VCF
    only applies to a small subset of cohorts
    """

    def expected_outputs(self, dataset: Dataset):
        """
        Expected to generate a VCF from the single-dataset MT
        """
        return {
            'vcf': (
                dataset.prefix()
                / 'vcf'
                / f'{get_workflow().output_version}-{dataset.name}.vcf.bgz'
            ),
            'index': (
                dataset.prefix()
                / 'vcf'
                / f'{get_workflow().output_version}-{dataset.name}.vcf.bgz.tbi'
            ),
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        only run this on manually defined list of cohorts
        """

        # only run this selectively, most datasets it's not required
        eligible_datasets = get_config()['workflow']['write_vcf']
        if dataset.name not in eligible_datasets:
            return None

        mt_path = inputs.as_path(target=dataset, stage=AnnotateDataset, key='mt')

        job = cohort_to_vcf_job(
            b=get_batch(),
            mt_path=mt_path,
            out_vcf_path=self.expected_outputs(dataset)['vcf'],
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset),
        )

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=job)


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
    analysis_keys=['index_name'],
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

        script = (
            f'cpg_workflows/dataproc_scripts/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {index_name} '
            f'--done-flag-path {done_flag_path} '
            f'--es-password {es_password()}'
        )
        pyfiles = ['seqr-loading-pipelines/hail_scripts']
        job_name = f'{dataset.name}: create ES index'

        if cluster_id := get_config()['hail'].get('dataproc', {}).get('cluster_id'):
            # noinspection PyProtectedMember
            j = dataproc._add_submit_job(
                batch=get_batch(),
                cluster_id=cluster_id,
                script=script,
                pyfiles=pyfiles,
                job_name=job_name,
                region='australia-southeast1',
            )
        else:
            j = dataproc.hail_dataproc_job(
                get_batch(),
                script,
                max_age='48h',
                packages=[
                    'cpg_workflows',
                    'elasticsearch==8.*',
                    'google',
                    'fsspec',
                    'gcloud',
                ],
                num_workers=2,
                num_secondary_workers=0,
                job_name=job_name,
                depends_on=inputs.get_jobs(dataset),
                scopes=['cloud-platform'],
                pyfiles=pyfiles,
                init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
            )
        j._preemptible = False
        j.attributes = (j.attributes or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]
        return self.make_outputs(dataset, data=index_name, jobs=jobs)
