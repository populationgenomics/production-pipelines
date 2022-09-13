#!/usr/bin/env python3

"""
Hail Query stages for the Seqr loader workflow.
"""

from cpg_utils import to_path
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
)

from jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs

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
            out_mt_path=self.expected_outputs(cohort)['mt'],
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
