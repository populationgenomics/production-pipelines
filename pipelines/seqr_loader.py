#!/usr/bin/env python3

"""
Batch pipeline to load data into seqr.
"""

import logging
import time

import click
import pandas as pd
from analysis_runner import dataproc
from cpg_pipes.storage import Path

from cpg_pipes import buckets, ref_data, utils
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.dataset import Cohort
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.pipeline import stage, Pipeline, CohortStage, DatasetStage
from cpg_pipes.pipeline.stage import StageInput, StageOutput
from cpg_pipes.stages.joint_genotyping import JointGenotypingStage
from cpg_pipes.stages.vqsr import VqsrStage
from cpg_pipes.storage import Namespace

logger = logging.getLogger(__file__)


def get_anno_tmp_bucket(cohort: Cohort) -> Path:
    """
    Path to write Hail Query intermediate files
    """
    return cohort.analysis_dataset.get_tmp_bucket() / 'mt'


@stage(required_stages=[JointGenotypingStage, VqsrStage])
class AnnotateCohortStage(CohortStage):
    """
    Re-annotate the entire cohort, including datasets that are not going to be loaded 
    """
    def expected_result(self, cohort: Cohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return get_anno_tmp_bucket(cohort) / 'combined.mt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        checkpoints_bucket = get_anno_tmp_bucket(cohort) / 'checkpoints'

        vcf_path = inputs.as_path(target=cohort, stage=JointGenotypingStage, id='vcf')
        annotated_siteonly_vcf_path = inputs.as_path(target=cohort, stage=VqsrStage)

        expected_path = self.expected_result(cohort)
        j = dataproc.hail_dataproc_job(
            self.b,
            f'{utils.QUERY_SCRIPTS_DIR}/seqr/vcf_to_mt.py '
            f'--vcf-path {vcf_path} '
            f'--site-only-vqsr-vcf-path {annotated_siteonly_vcf_path} '
            f'--dest-mt-path {expected_path} '
            f'--bucket {checkpoints_bucket} '
            f'--disable-validation '
            f'--make-checkpoints '
            + ('--overwrite ' if not self.check_intermediates else ''),
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            job_name='Make MT and annotate cohort',
            depends_on=inputs.get_jobs(),
            # Default Hail's VEP initialization script (triggered by --vep) 
            # installs VEP=v95; if we want v105, we have to use a modified 
            # vep-GRCh38.sh (with --init) from production-pipelines/vep/vep-GRCh38.sh
            init=['gs://cpg-reference/vep/vep-GRCh38.sh'],
            worker_machine_type='n1-highmem-8',
            worker_boot_disk_size=200,
            secondary_worker_boot_disk_size=200,
            num_secondary_workers=50,
            num_workers=8,
        )
        return self.make_outputs(cohort, data=expected_path, jobs=[j])


@stage(required_stages=[AnnotateCohortStage])
class AnnotateDatasetStage(DatasetStage):
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """
    def expected_result(self, dataset: Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        return dataset.get_analysis_bucket() / 'mt' / f'{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        annotated_mt_path = inputs.as_path(
            target=dataset.cohort, 
            stage=AnnotateCohortStage
        )

        # Make a list of dataset samples to subset from the entire matrix table
        sample_ids = [s.id for s in dataset.get_samples()]
        proj_tmp_bucket = dataset.get_tmp_bucket()
        subset_path = proj_tmp_bucket / 'seqr-samples.txt'
        with subset_path.open('w') as f:
            f.write('\n'.join(sample_ids))

        expected_path = self.expected_result(dataset)    
        j = dataproc.hail_dataproc_job(
            self.b,
            f'{utils.QUERY_SCRIPTS_DIR}/seqr/subset_mt.py '
            f'--mt-path {annotated_mt_path} '
            f'--out-mt-path {expected_path} '
            f'--subset-tsv {subset_path}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=20,
            num_workers=5,
            job_name=f'{dataset.name}: annotate dataset',
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(dataset, data=expected_path, jobs=[j])


@stage(required_stages=[AnnotateDatasetStage])
class LoadToEsStage(DatasetStage):
    """
    Create a Seqr index.
    """
    def expected_result(self, dataset: Dataset) -> None:
        """
        Expected to generate a Seqr index, which is not a file
        """
        return None

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetStage)
        version = time.strftime('%Y%m%d-%H%M%S')

        j = dataproc.hail_dataproc_job(
            self.b,
            f'{utils.QUERY_SCRIPTS_DIR}/seqr/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {dataset.name}-{version} '
            f'--es-index-min-num-shards 1 '
            f'{"--prod" if dataset.namespace == Namespace.MAIN else ""}',
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=2,
            job_name=f'{dataset.name}: create ES index',
            depends_on=inputs.get_jobs(),
            scopes=['cloud-platform'],
        )
        return self.make_outputs(dataset, jobs=[j])


@stage(required_stages=[AnnotateDatasetStage])
class ToVcfStage(DatasetStage):
    """
    Convers the annotated matrix table to a pVCF
    """
    def expected_result(self, dataset: Dataset) -> Path:
        """
        Generates an indexed VCF
        """
        return dataset.get_analysis_bucket() / 'vcf' / f'{dataset.name}.vcf.bgz'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetStage)

        j = dataproc.hail_dataproc_job(
            self.b,
            f'{utils.QUERY_SCRIPTS_DIR}/seqr/mt_to_vcf.py '
            f'--mt-path {dataset_mt_path} '
            f'--out-vcf-path {self.expected_result(dataset)}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=20,
            num_workers=5,
            job_name=f'{dataset.name}: export to VCF',
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(dataset, self.expected_result(dataset), jobs=[j])


def _make_seqr_metadata_files(
    dataset: Dataset, 
    bucket: Path,
    local_dir: Path, 
    overwrite: bool = False
):
    """
    Create Seqr metadata files
    """
    samplemap_bucket_path = bucket / 'seqr' / f'{dataset.name}-sample-map.csv'
    igv_paths_path = local_dir / f'{dataset.name}-igv-paths.tsv'

    # Sample map
    if not buckets.can_reuse(samplemap_bucket_path, overwrite):
        df = pd.DataFrame({
            'cpg_id': s.id,
            'individual_id': s.participant_id,
        } for s in dataset.get_samples())
        df.to_csv(str(samplemap_bucket_path), sep=',', index=False, header=False)

    # IGV
    if not buckets.can_reuse(igv_paths_path, overwrite):
        df = pd.DataFrame({
            'individual_id': s.participant_id,
            'cram_path': s.get_cram_path(),
            'cram_sample_id': s.id,
        } for s in dataset.get_samples() if s.get_cram_path())
        df.to_csv(str(igv_paths_path), sep='\t', index=False, header=False)

    logger.info(f'Seqr sample map: {samplemap_bucket_path}')
    logger.info(f'IGV seqr paths: {igv_paths_path}')


@click.command()
@click.option(
    '--hc-shards-num',
    'hc_shards_num',
    type=click.INT,
    default=ref_data.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
    help='Number of intervals to devide the genome for gatk HaplotypeCaller',
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
    '--make-seqr-metadata/--no-make-seqr-metadata',
    'make_seqr_metadata',
    default=True,
    is_flag=True,
    help='Make Seqr metadata',
)
@click.option(
    '--ped-checks/--no-ped-checks',
    'ped_checks',
    default=True,
    is_flag=True,
    help='Perform fingerprinting and PED checks',
)
@pipeline_click_options
def main(
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    ped_checks: bool,
    **kwargs,
):
    """
    Entry point, decorated by pipeline click options.
    """
    make_seqr_metadata = kwargs.pop('make_seqr_metadata')

    pipeline = Pipeline(
        name='seqr_loader',
        description='Seqr loader',
        config=dict(
            ped_checks=ped_checks,
            hc_shards_num=hc_shards_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
        ),
        **kwargs,
    )
    pipeline.submit_batch()

    if make_seqr_metadata:
        for dataset in pipeline.cohort.get_datasets():
            _make_seqr_metadata_files(
                dataset=dataset,
                bucket=pipeline.cohort.analysis_dataset.get_analysis_bucket(
                    version=pipeline.version,
                ),
                local_dir=pipeline.local_dir,
            )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
