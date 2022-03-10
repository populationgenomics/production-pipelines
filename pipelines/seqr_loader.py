#!/usr/bin/env python3

"""
Batch pipeline to laod data into seqr
"""

import logging
import time
from os.path import join, basename
from typing import Optional, List

import click
import hail as hl
import pandas as pd

from analysis_runner import dataproc

from cpg_pipes import buckets, ref_data, utils
from cpg_pipes.jobs import align, split_intervals, haplotype_caller, \
    pedigree, fastqc
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.namespace import Namespace
from cpg_pipes.pipeline.analysis import AnalysisType, GvcfPath, CramPath
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.pipeline import stage, Pipeline, PipelineError
from cpg_pipes.pipeline.stage import SampleStage, CohortStage, DatasetStage, \
    StageInput, StageOutput
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.cohort import Cohort

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage(skipped=True)
class FastqcStage(SampleStage):
    """
    Run FastQC on alignment inputs
    """
    def expected_result(self, sample: Sample):
        """
        Stage is expected to generate a fastqc report
        """
        return f'{sample.dataset.get_bucket()}/qc/fastqc/{sample.id}.html'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "fastqc" function implemented in the jobs module
        """
        if not sample.alignment_input:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(
                    f'No alignment input found for {sample.id}. '
                    f'Checked: Sequence entry and type=CRAM Analysis entry'
                )

        job = fastqc.fastqc(
            b=self.pipe.b,
            output_fpath=self.expected_result(sample),
            alignment_input=sample.alignment_input,
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[job]
        )
    

@stage(sm_analysis_type=AnalysisType.CRAM)
class CramStage(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """
    def expected_result(self, sample: Sample):
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        return sample.cram_path.path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "align" function implemented in the jobs module
        """
        if not sample.alignment_input:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(
                    f'No alignment input found for {sample.id}. '
                    f'Checked: Sequence entry and type=CRAM Analysis entry'
                )

        cram_job = align.align(
            b=self.pipe.b,
            alignment_input=sample.alignment_input,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
            overwrite=not self.pipe.check_intermediates,
            smdb=self.pipe.get_db(),
            prev_batch_jobs=self.pipe.prev_batch_jobs,
            number_of_shards_for_realignment=(
                10 if isinstance(sample.alignment_input, CramPath) else None
            )
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[cram_job]
        )


@stage(required_stages=CramStage)
class CramSomalierStage(SampleStage):
    """
    Genereate fingerprints from CRAMs for pedigree checks
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate the fingerprints file
        """
        return sample.cram_path.somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using a function from the jobs module
        """
        cram_path = inputs.as_path(target=sample, stage=CramStage)
        expected_path = self.expected_result(sample)
        j, _ = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=CramPath(cram_path),
            out_fpath=expected_path,
            overwrite=not self.pipe.check_intermediates,
            label='(CRAMs)',
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage(required_stages=CramSomalierStage)
class CramPedCheckStage(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints
    """
    def expected_result(self, dataset: Dataset):
        """
        We don't expect any output - just rerun it every time
        """

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = inputs.as_path_by_target(stage=CramSomalierStage)

        j, somalier_samples_path, _ = pedigree.add_pedigree_jobs(
            self.pipe.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(CRAMs)',
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=somalier_samples_path, jobs=[j])


@stage(required_stages=CramStage, sm_analysis_type=AnalysisType.GVCF)
class GvcfStage(SampleStage):
    """
    Use HaplotypeCaller to genotype individual samples
    """
    hc_intervals = None

    def expected_result(self, sample: Sample):
        """
        Generate a GVCF and corresponding TBI index
        """
        return sample.gvcf_path.path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Use function from the jobs module
        """
        hc_shards_num = self.pipe.config.get('hc_shards_num', 1)
        if GvcfStage.hc_intervals is None and hc_shards_num > 1:
            GvcfStage.hc_intervals = split_intervals.get_intervals(
                b=self.pipe.b,
                scatter_count=hc_shards_num,
            )
        gvcf_job = haplotype_caller.produce_gvcf(
            b=self.pipe.b,
            output_path=self.expected_result(sample),
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
            cram_path=sample.cram_path,
            intervals=GvcfStage.hc_intervals,
            number_of_intervals=hc_shards_num,
            tmp_bucket=self.pipe.tmp_bucket,
            overwrite=not self.pipe.check_intermediates,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[gvcf_job]
        )


@stage(required_stages=GvcfStage)
class GvcfSomalierStage(SampleStage):
    """
    Genereate fingerprints from GVCFs for pedigree checks
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate the fingerprints file
        """
        return sample.gvcf_path.somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Use function from pedigree module
        """
        gvcf_path = inputs.as_path(target=sample, stage=GvcfStage)
        expected_path = self.expected_result(sample)
        j, _ = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=GvcfPath(gvcf_path),
            out_fpath=expected_path,
            overwrite=not self.pipe.check_intermediates,
            label='(GVCFs)',
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage(required_stages=GvcfSomalierStage)
class GvcfPedCheckStage(DatasetStage):
    """
    Check pedigree from GVCFs
    """
    def expected_result(self, dataset: Dataset):
        """
        We don't expect any output - just rerun it every time
        """

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Use function from pedigree module
        """
        fp_by_sid = inputs.as_path_by_target(stage=GvcfSomalierStage)

        j, somalier_samples_path, _ = pedigree.add_pedigree_jobs(
            self.pipe.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(GVCFs)',
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=somalier_samples_path, jobs=[j])


@stage(required_stages=GvcfStage, sm_analysis_type=AnalysisType.JOINT_CALLING)
class JointGenotypingStage(CohortStage):
    """
    Joint-calling of GVCFs together
    """
    def expected_result(self, cohort: Cohort):
        """
        Generate a pVCF and a site-only VCF
        """
        samples_hash = utils.hash_sample_ids(cohort.get_all_sample_ids())
        expected_jc_vcf_path = (
            f'{self.pipe.tmp_bucket}/joint_calling/{samples_hash}.vcf.gz'
        )
        return {
            'vcf': expected_jc_vcf_path,
            'siteonly': expected_jc_vcf_path.replace('.vcf.gz', '-siteonly.vcf.gz'),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Use function defined in jobs module
        """
        gvcf_by_sid = {
            k: GvcfPath(v) for k, v 
            in inputs.as_path_by_target(stage=GvcfStage).items()
        }

        not_found_gvcfs: List[str] = []
        for sid, gvcf_path in gvcf_by_sid.items():
            if gvcf_path is None:
                logger.error(f'Joint genotyping: could not find GVCF for {sid}')
                not_found_gvcfs.append(sid)
        if not_found_gvcfs:
            raise PipelineError(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting'
            )

        jc_job = make_joint_genotyping_jobs(
            b=self.pipe.b,
            out_vcf_path=self.expected_result(cohort)['vcf'],
            out_siteonly_vcf_path=self.expected_result(cohort)['siteonly'],
            samples=cohort.get_all_samples(),
            genomicsdb_bucket=f'{self.pipe.analysis_bucket}/genomicsdbs',
            tmp_bucket=self.pipe.tmp_bucket,
            gvcf_by_sid=gvcf_by_sid,
            overwrite=not self.pipe.check_intermediates,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
            tool=JointGenotyperTool.GnarlyGenotyper 
            if self.pipe.config.get('use_gnarly', False) 
            else JointGenotyperTool.GenotypeGVCFs,
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(
            cohort, 
            data=self.expected_result(cohort), 
            jobs=[jc_job]
        )


@stage(required_stages=JointGenotypingStage)
class VqsrStage(CohortStage):
    """
    Variant filtering of joint-called VCF
    """
    def expected_result(self, cohort: Cohort):
        """
        Expects to generate one site-only VCF
        """
        samples_hash = utils.hash_sample_ids(cohort.get_all_sample_ids())
        expected_jc_vcf_path = (
            f'{self.pipe.tmp_bucket}/vqsr/{samples_hash}-siteonly.vcf.gz'
        )
        return expected_jc_vcf_path

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Use function defined in jobs module
        """
        siteonly_vcf_path = inputs.as_path(
            stage=JointGenotypingStage, target=cohort, id='siteonly'
        )

        tmp_vqsr_bucket = f'{self.pipe.tmp_bucket}/vqsr'
        logger.info(f'Queueing VQSR job')
        expected_path = self.expected_result(cohort)
        vqsr_job = make_vqsr_jobs(
            b=self.pipe.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            work_bucket=tmp_vqsr_bucket,
            gvcf_count=len(cohort.get_all_samples()),
            depends_on=inputs.get_jobs(),
            output_vcf_path=expected_path,
            use_as_annotations=self.pipe.config.get('use_as_vqsr', True),
            overwrite=not self.pipe.check_intermediates,
        )
        return self.make_outputs(cohort, data=expected_path, jobs=[vqsr_job])


def get_anno_tmp_bucket(pipe: Pipeline):
    """
    Path to write Hail Query intermediate files
    """
    return join(pipe.tmp_bucket, 'mt')


@stage(required_stages=[JointGenotypingStage, VqsrStage])
class AnnotateCohortStage(CohortStage):
    """
    Re-annotate the entire cohort, including datasets that are not going to be loaded 
    """
    def expected_result(self, cohort: Cohort):
        """
        Expected to write a matrix table
        """
        return join(get_anno_tmp_bucket(self.pipe), 'combined.mt')

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        checkpoints_bucket = join(get_anno_tmp_bucket(self.pipe), 'checkpoints')

        vcf_path = inputs.as_path(target=cohort, stage=JointGenotypingStage, id='vcf')
        annotated_siteonly_vcf_path = inputs.as_path(target=cohort, stage=VqsrStage)

        expected_path = self.expected_result(cohort)
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "vcf_to_mt.py")} '
            f'--vcf-path {vcf_path} '
            f'--site-only-vqsr-vcf-path {annotated_siteonly_vcf_path} '
            f'--dest-mt-path {expected_path} '
            f'--bucket {checkpoints_bucket} '
            f'--disable-validation '
            f'--make-checkpoints '
            + ('--overwrite ' if not self.pipe.check_intermediates else ''),
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
    def expected_result(self, dataset: Dataset):
        """
        Expected to generate a matrix table
        """
        return f'{self.pipe.analysis_bucket}/mt/{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        output_datasets = self.pipe.config.get('output_datasets', self.pipe.get_datasets())
        if dataset.stack not in output_datasets:
            logger.info(
                f'Skipping annotating dataset {dataset.stack} because it is not'
                f'in the --output-datasets: {output_datasets}'
            )
            return self.make_outputs(dataset)
        
        annotated_mt_path = inputs.as_path(
            target=self.pipe.cohort, 
            stage=AnnotateCohortStage
        )

        # Make a list of dataset samples to subset from the entire matrix table
        sample_ids = [s.id for s in dataset.get_samples()]
        proj_tmp_bucket = dataset.get_tmp_bucket()
        subset_path = f'{proj_tmp_bucket}/seqr-samples.txt'
        with hl.hadoop_open(subset_path, 'w') as f:
            f.write('\n'.join(sample_ids))

        expected_path = self.expected_result(dataset)    
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "subset_mt.py")} '
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
    Create a Seqr index for requested "output_datasets"
    """
    def expected_result(self, dataset: Dataset):
        """
        Expected to generate a Seqr index, which is not a file
        """
        return None

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        output_datasets = self.pipe.config.get(
            'output_datasets', 
            self.pipe.get_datasets()
        )
        if dataset.stack not in output_datasets:
            logger.info(
                f'Skipping loading a dataset {dataset.stack} because it is not'
                f'in the --output-datasets: {output_datasets}'
            )
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetStage)

        timestamp = time.strftime('%Y%m%d-%H%M%S')
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_es.py")} '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {dataset.name}-{self.pipe.output_version}-{timestamp} '
            f'--es-index-min-num-shards 1 '
            f'{"--prod" if self.pipe.namespace == Namespace.MAIN else ""}',
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
    def expected_result(self, dataset: Dataset):
        """
        Generates an indexed VCF
        """
        return f'{self.pipe.analysis_bucket}/mt/{dataset.name}.vcf.bgz'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        output_datasets = self.pipe.config.get('output_datasets', self.pipe.get_datasets())
        if dataset.stack not in output_datasets:
            logger.info(
                f'Skipping loading dataset {dataset.stack} because it is not'
                f'in the --output-datasets: {output_datasets}'
            )
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetStage)

        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_vcf.py")} '
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


def make_pipeline(
    input_datasets: List[str],
    output_datasets: Optional[List[str]],
    output_version: str,
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    **kwargs,
) -> Pipeline:
    """
    Create the seqr-loader Pipeline
    """
    assert input_datasets
    title = f'Seqr loading: joint call from: {", ".join(input_datasets)}'
    if output_datasets:
        title += f', ES index for: {", ".join(output_datasets)}'
    title += f', version {output_version}'

    if not output_datasets:
        output_datasets = input_datasets
    if not all(op in input_datasets for op in output_datasets):
        raise click.BadParameter(
            f'All output datasets must be contained within the specified input '
            f'datasets. Input dataset: {input_datasets}, output datasets: '
            f'{output_datasets}'
        )

    pipeline = Pipeline(
        name='seqr_loader',
        description=title,
        config=dict(
            output_datasets=output_datasets,
            skip_ped_checks=skip_ped_checks,
            hc_shards_num=hc_shards_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
        ),
        input_datasets=input_datasets,
        output_version=output_version,
        **kwargs,
    )
    return pipeline


def _make_seqr_metadata_files(
    dataset: Dataset, 
    bucket: str, 
    local_dir: str, 
    overwrite: bool = False
):
    """
    Create Seqr metadata files
    """
    sample_map_path = join(local_dir, f'{dataset.name}-sample-map.csv')
    igv_paths_path = join(local_dir, f'{dataset.name}-igv-paths.tsv')

    # Sample map
    if not buckets.can_reuse(sample_map_path, overwrite):
        df = pd.DataFrame({
            'cpg_id': s.id,
            'individual_id': s.participant_id,
        } for s in dataset.get_samples())
        df.to_csv(sample_map_path, sep=',', index=False, header=False)
        # Sample map has to sit on a bucket
        sample_map_path = buckets.gsutil_cp(
            sample_map_path, 
            f'{bucket}/seqr/{basename(sample_map_path)}'
        )

    # IGV
    if not buckets.can_reuse(igv_paths_path, overwrite):
        df = pd.DataFrame({
            'individual_id': s.participant_id,
            'cram_path': s.cram_path,
            'cram_sample_id': s.id,
        } for s in dataset.get_samples() if s.cram_path)
        df.to_csv(igv_paths_path, sep='\t', index=False, header=False)

    logger.info(f'Seqr sample map: {sample_map_path}')
    logger.info(f'IGV seqr paths: {igv_paths_path}')


@click.command()
@pipeline_click_options
@click.option(
    '--output-dataset',
    'output_datasets',
    multiple=True,
    help='Only create ES indicies for the dataset(s). Can be set multiple times. '
    'Defaults to --input-datasets. The name of the ES index will be suffixed '
    'with the dataset version (set by --version)',
)
@click.option(
    '--skip-ped-checks',
    'skip_ped_checks',
    is_flag=True,
    help='Skip checking provided sex and pedigree against the inferred one',
)
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
def main(**kwargs):  # pylint: disable=missing-function-docstring
    make_seqr_metadata = kwargs.pop('make_seqr_metadata')
    
    pipeline = make_pipeline(**kwargs)
    pipeline.submit_batch()

    if make_seqr_metadata:
        for dataset in pipeline.cohort.get_datasets():
            if dataset.stack in pipeline.config.get('output_datasets', []):
                _make_seqr_metadata_files(
                    dataset=dataset,
                    bucket=pipeline.analysis_bucket,
                    local_dir=pipeline.local_dir,
                )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
