#!/usr/bin/env python3

"""
Hail Query stages for the Seqr loader workflow.
"""
from typing import Any

from google.api_core.exceptions import PermissionDenied

from cpg_utils import Path
from cpg_utils.cloud import read_secret
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.jobs.seqr_loader import annotate_dataset_jobs, cohort_to_vcf_job
from cpg_workflows.query_modules import seqr_loader
from cpg_workflows.targets import Dataset, MultiCohort
from cpg_workflows.utils import get_logger, tshirt_mt_sizing
from cpg_workflows.workflow import (
    DatasetStage,
    MultiCohortStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)

from .joint_genotyping import JointGenotyping
from .vep import Vep
from .vqsr import Vqsr


@stage(required_stages=[JointGenotyping, Vqsr, Vep])
class AnnotateCohort(MultiCohortStage):
    """
    Re-annotate the entire cohort.
    """

    def expected_outputs(self, multicohort: MultiCohort):
        """
        Expected to write a matrix table.
        """
        return {'mt': self.prefix / 'cohort.mt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Apply VEP and VQSR annotations to all-sample callset
        """
        outputs = self.expected_outputs(multicohort)
        vcf_path = inputs.as_path(target=multicohort, stage=JointGenotyping, key='vcf')
        siteonly_vqsr_vcf_path = inputs.as_path(target=multicohort, stage=Vqsr, key='siteonly')
        vep_ht_path = inputs.as_path(target=multicohort, stage=Vep, key='ht')

        j = get_batch().new_job('annotate cohort', self.get_job_attrs(multicohort))
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                seqr_loader,
                seqr_loader.annotate_cohort.__name__,
                str(vcf_path),
                str(outputs['mt']),
                str(vep_ht_path),
                str(siteonly_vqsr_vcf_path) if siteonly_vqsr_vcf_path else None,
                str(self.tmp_prefix / 'checkpoints'),
                setup_gcp=True,
            ),
        )
        if depends_on := inputs.get_jobs(multicohort):
            j.depends_on(*depends_on)

        return self.make_outputs(multicohort, data=outputs, jobs=j)


def _snv_es_index_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to es-index analysis object
    https://github.com/populationgenomics/metamist/issues/539
    """
    return {'seqr-dataset-type': 'VARIANTS'}


@stage(required_stages=[AnnotateCohort], analysis_type='custom', analysis_keys=['mt'])
class AnnotateDataset(DatasetStage):
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """

    def expected_outputs(self, dataset: Dataset):
        """
        Expected to generate a matrix table
        """
        return {'mt': (dataset.prefix() / 'mt' / f'{get_workflow().output_version}-{dataset.name}.mt')}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Annotate MT with genotype and consequence data for Seqr
        """
        # only create dataset MTs for datasets specified in the config
        eligible_datasets = config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
        if dataset.name not in eligible_datasets:
            get_logger().info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return None

        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohort, key='mt')

        jobs = annotate_dataset_jobs(
            mt_path=mt_path,
            sequencing_group_ids=dataset.get_sequencing_group_ids(),
            out_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_prefix=self.tmp_prefix / dataset.name / 'checkpoints',
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset),
        )

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)


@stage(required_stages=[AnnotateDataset], analysis_type='custom', analysis_keys=['vcf'])
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
            'vcf': (dataset.prefix() / 'vcf' / f'{get_workflow().output_version}-{dataset.name}.vcf.bgz'),
            'index': (dataset.prefix() / 'vcf' / f'{get_workflow().output_version}-{dataset.name}.vcf.bgz.tbi'),
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Run a MT -> VCF extraction on selected cohorts
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
    update_analysis_meta=_snv_es_index_meta,
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
        index_name = f'{dataset.name}-{sequencing_type}-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Transforms the MT into a Seqr index, no DataProc
        """
        # only create the elasticsearch index for the datasets specified in the config
        eligible_datasets = config_retrieve(['workflow', 'create_es_index_for_datasets'], default=[])
        if dataset.name not in eligible_datasets:
            get_logger().info(f'Skipping ES index creation for {dataset}')
            return None

        # try to generate a password here - we'll find out inside the script anyway, but
        # by that point we'd already have localised the MT, wasting time and money
        try:
            _es_password_string = es_password()
        except PermissionDenied:
            get_logger().warning(f'No permission to access ES password, skipping for {dataset}')
            return self.make_outputs(dataset)
        except KeyError:
            get_logger().warning(f'ES section not in config, skipping for {dataset}')
            return self.make_outputs(dataset)

        # get the absolute path to the MT
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDataset, key='mt'))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        outputs = self.expected_outputs(dataset)

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_bash_job(f'Generate {index_name} from {mt_path}')
        if config_retrieve(['workflow', 'es_index', 'spot_instance'], default=True) is False:
            # Use a non-preemptible instance if spot_instance is False in the config
            job = job.spot(is_spot=False)

        req_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )

        job.cpu(4).storage(f'{req_storage}Gi').memory('lowmem').image(config_retrieve(['workflow', 'driver_image']))

        # localise the MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

        # run the export from the localised MT - this job writes no new data, just transforms and exports over network
        job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
