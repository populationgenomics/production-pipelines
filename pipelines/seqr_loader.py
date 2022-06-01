#!/usr/bin/env python3

"""
Seqr loading pipeline: FASTQ -> ElasticSearch index.
"""

import logging
import os
import time

import click
import yaml

from cpg_pipes.pipeline.cli_opts import choice_from_enum, val_to_enum, create_pipeline
from cpg_pipes.types import SequencingType
from google.cloud import secretmanager

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from cpg_pipes.pipeline import (
    pipeline_options,
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
)
from cpg_pipes.stages.vep import Vep
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.targets import Cohort, Dataset

logger = logging.getLogger(__file__)

SUPPORTED_SEQUENCING_TYPES = [
    SequencingType.GENOME,
]


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
            hail_billing_project=self.hail_billing_project,
            hail_bucket=self.hail_bucket,
            overwrite=not self.check_intermediates,
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
            hail_billing_project=self.hail_billing_project,
            hail_bucket=self.hail_bucket,
            job_attrs=self.get_job_attrs(dataset),
            overwrite=not self.check_intermediates,
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
            es_datasets := self.pipeline_config.get('create_es_index_for_datasets')
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


@click.command()
@click.option(
    '--es-dataset',
    'create_es_index_for_datasets',
    multiple=True,
    help=f'Create Seqr ElasticSearch indices for these datasets.',
)
@click.option(
    '--realign',
    '--realign-from-cram-version',
    'realign_from_cram_version',
    help='Realign CRAM whenever available, instead of using FASTQ. '
         'The parameter value should correspond to CRAM version '
         '(e.g. v0 in gs://cpg-fewgenomes-main/cram/v0/CPG0123.cram)'
)
@click.option(
    '--use-gnarly/--no-use-gnarly',
    'use_gnarly',
    default=False,
    is_flag=True,
    help='Use GnarlyGenotyper instead of GenotypeGVCFs',
)
@click.option(
    '--use-as-vqsr/--no-use-as-vqsr',
    'use_as_vqsr',
    default=True,
    is_flag=True,
    help='Use allele-specific annotations for VQSR',
)
@click.option(
    '--cram-qc/--no-cram-qc',
    'cram_qc',
    default=True,
    is_flag=True,
    help='Run CRAM QC and PED checks',
)
@click.option(
    '--exome-bed',
    'exome_bed',
    help=f'BED file with exome regions',
)
@click.option(
    '--sequencing-type',
    'sequencing_type',
    type=choice_from_enum(SequencingType),
    callback=val_to_enum(SequencingType),
    help='Limit to data with this sequencing type',
    default=SequencingType.GENOME.value,
)
@pipeline_options
def main(
    datasets: list[str],
    create_es_index_for_datasets: list[str],
    sequencing_type: SequencingType,
    **kwargs,
):
    """
    Seqr loading pipeline: FASTQ -> ElasticSearch index.
    """
    if not datasets:
        # Parsing dataset names from the analysis-runner Seqr stack:
        from urllib import request

        seqr_stack_url = (
            'https://raw.githubusercontent.com/populationgenomics/analysis-runner/main'
            '/stack/Pulumi.seqr.yaml'
        )
        with request.urlopen(seqr_stack_url) as f:
            value = yaml.safe_load(f)['config']['datasets:depends_on']
            datasets = [d.strip('"') for d in value.strip('[] ').split(', ')]
            logger.info(f'Found Seqr datasets: {datasets}')

    description = 'Seqr Loader'
    if v := kwargs.get('version'):
        description += f', {v}'
    active_datasets = set(datasets) - set(kwargs.get('skip_datasets', []))
    description += f': [{", ".join(sorted(active_datasets))}]'
    if create_es_index_for_datasets:
        description += f' → [{", ".join(sorted(create_es_index_for_datasets))}]'

    if sequencing_type not in SUPPORTED_SEQUENCING_TYPES:
        raise click.BadParameter(
            f'Unsupported sequencing data type {sequencing_type.value}. '
            f'Supported types: {[st.value for st in SUPPORTED_SEQUENCING_TYPES]} '
        )

    kwargs['analysis_dataset'] = kwargs.get('analysis_dataset', 'seqr')
    kwargs['name'] = kwargs.get('name') or 'Seqr Loader'
    kwargs['description'] = description
    kwargs['datasets'] = kwargs.get('datasets') or datasets
    kwargs['create_es_index_for_datasets'] = create_es_index_for_datasets
    kwargs['sequencing_type'] = sequencing_type
    pipeline = create_pipeline(**kwargs)
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
