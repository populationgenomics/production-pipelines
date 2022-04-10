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
from analysis_runner import dataproc

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from cpg_pipes.pipeline import (
    pipeline_click_options,
    stage,
    create_pipeline,
    StageInput,
    StageOutput,
    CohortStage,
    DatasetStage,
)
from cpg_pipes.refdata import RefData
from cpg_pipes.stages.vep import VepStage
from cpg_pipes.stages.joint_genotyping import JointGenotypingStage
from cpg_pipes.stages.vqsr import VqsrStage
from cpg_pipes.targets import Cohort, Dataset

logger = logging.getLogger(__file__)


@stage(required_stages=[JointGenotypingStage, VepStage, VqsrStage])
class AnnotateCohort(CohortStage):
    """
    Re-annotate the entire cohort.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Expected to write a matrix table.
        """
        h = cohort.alignment_inputs_hash()
        return cohort.analysis_dataset.get_analysis_bucket() / 'mt' / f'{h}.mt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotypingStage, id='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(target=cohort, stage=VqsrStage)
        vep_ht_path = inputs.as_path(target=cohort, stage=VepStage)

        mt_path = self.expected_outputs(cohort)

        jobs = annotate_cohort_jobs(
            b=self.b,
            vcf_path=vcf_path,
            vep_ht_path=vep_ht_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_vcf_path,
            output_mt_path=mt_path,
            checkpoints_bucket=self.tmp_bucket / 'seqr_loader' / 'checkpoints',
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
        h = dataset.cohort.alignment_inputs_hash()
        return self.tmp_bucket / f'{h}-{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        mt_path = inputs.as_path(target=dataset.cohort, stage=AnnotateCohort)

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
class LoadToEsStage(DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: Dataset) -> None:
        """
        Expected to generate a Seqr index, which is not a file
        """
        return None

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        if (
            'output_datasets' in self.pipeline_config
            and dataset.name not in self.pipeline_config['output_datasets']
        ):
            # Skipping dataset that wasn't explicitly requested to upload to ES:
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDataset)
        version = time.strftime('%Y%m%d-%H%M%S')

        j = dataproc.hail_dataproc_job(
            self.b,
            f'{utils.QUERY_SCRIPTS_DIR}/seqr/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-host elasticsearch.es.australia-southeast1.gcp.elastic-cloud.com '
            f'--es-port 9243 '
            f'--es-username seqr '
            f'--es-password {_read_es_password()} '
            f'--es-index {dataset.name}-{version} '
            f'--es-index-min-num-shards 1 '
            f'--use-spark ',  # es export doesn't work with the service backend
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=2,
            job_name=f'{dataset.name}: create ES index',
            depends_on=inputs.get_jobs(dataset),
            scopes=['cloud-platform'],
        )
        j.attributes = (self.get_job_attrs(dataset),)
        return self.make_outputs(dataset, jobs=[j])


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
    '--output-dataset',
    'output_datasets',
    multiple=True,
    help=f'Datasets to load into Seqr',
)
@click.option(
    '--hc-intervals-num',
    'hc_intervals_num',
    type=click.INT,
    default=RefData.number_of_haplotype_caller_intervals,
    help='Number of intervals to devide the genome for sample genotyping with '
    'gatk HaplotypeCaller',
)
@click.option(
    '--jc-intervals-num',
    'jc_intervals_num',
    type=click.INT,
    default=RefData.number_of_joint_calling_intervals,
    help='Number of intervals to devide the genome for joint genotyping with GATK',
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
    '--ped-checks/--no-ped-checks',
    'ped_checks',
    default=True,
    is_flag=True,
    help='Perform fingerprinting and PED checks',
)
@click.option(
    '--exome-bed',
    'exome_bed',
    help=f'BED file with exome regions',
)
@pipeline_click_options
def main(
    datasets: list[str],
    output_datasets: list[str],
    hc_intervals_num: int,
    jc_intervals_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    ped_checks: bool,
    exome_bed: str,
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

    description = 'Seqr Loader'
    if v := kwargs.get('version'):
        description += f', {v}'
    input_datasets = set(datasets) - set(kwargs.get('skip_datasets', []))
    description += f': [{", ".join(input_datasets)}]'
    if output_datasets:
        description += f' → [{", ".join(output_datasets)}]'

    kwargs['analysis_dataset'] = 'seqr'
    pipeline = create_pipeline(
        name='seqr_loader',
        description=description,
        datasets=datasets,
        config=dict(
            ped_checks=ped_checks,
            hc_intervals_num=hc_intervals_num,
            jc_intervals_num=jc_intervals_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
            exome_bed=exome_bed,
            output_datasets=output_datasets,
        ),
        **kwargs,
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
