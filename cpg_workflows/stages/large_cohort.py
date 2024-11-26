import logging
from functools import cache
from itertools import chain
from logging import config
from typing import TYPE_CHECKING, Any, Final, Tuple

from sympy import root

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, genome_build, get_config, image_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.jobs.gcta_PCA import run_PCA
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.utils import slugify, tshirt_mt_sizing
from cpg_workflows.workflow import (
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

if TYPE_CHECKING:
    from graphql import DocumentNode

    from hailtop.batch.job import BashJob, PythonJob


HAIL_QUERY: Final = 'hail query'


@cache
def relatedness_version() -> str:
    """
    generate the relatedness version
    """
    if not (relatedness_str := config_retrieve(['large_cohort', 'output_versions', 'relatedness'], default=None)):
        return get_workflow().output_version
    return slugify(relatedness_str)


@cache
def sample_qc_version() -> str:
    """
    generate the sample_qc version
    """
    if not (sample_qc_str := config_retrieve(['large_cohort', 'output_versions', 'sample_qc'], default=None)):
        return get_workflow().output_version
    return slugify(sample_qc_str)


@cache
def vds_version() -> str:
    if not (vds_version_str := get_config()['workflow'].get('vds_version')):
        return get_workflow().output_version
    vds_version_str = slugify(vds_version_str)
    # if not vds_version_str.startswith('v'):
    #     vds_version_str = f'v{vds_version_str}'
    return vds_version_str


@cache
def dense_subset_version() -> str:
    """
    generate the dense_subset version
    """
    if not (dense_subset_str := config_retrieve(['large_cohort', 'output_versions', 'dense_subset'], default=None)):
        return get_workflow().output_version
    return slugify(dense_subset_str)


@cache
def dense_background_vds_version(dataset: str) -> str:
    """
    generate the dense_background version
    """
    dataset = dataset.replace('-', '_')
    if not (
        dense_bg_vds_version_str := config_retrieve(
            ['large_cohort', 'output_versions', f'{dataset}_dense_background'],
            default=None,
        )
    ):
        return get_workflow().output_version
    dense_bg_vds_version_str = slugify(dense_bg_vds_version_str)
    if not dense_bg_vds_version_str.startswith('v'):
        dense_bg_vds_version_str = f'v{dense_bg_vds_version_str}'
    return dense_bg_vds_version_str


@cache
def gcta_version() -> str:
    """
    generate the GCTA version
    """
    if not (dense_subset_str := config_retrieve(['large_cohort', 'output_versions', 'gcta'], default=None)):
        return get_workflow().output_version
    return slugify(dense_subset_str)


@stage(analysis_type='combiner')
class Combiner(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        workflow_config = config_retrieve('workflow')
        combiner_config = config_retrieve('combiner')
        output_vds_name: str = slugify(
            f"{workflow_config['cohort']}-{workflow_config['sequencing_type']}-{combiner_config['vds_version']}",
        )
        return cohort.analysis_dataset.prefix() / 'vds' / f'{output_vds_name}.vds'

    def get_vds_ids_output(self, vds_id: int) -> Tuple[str, list[str]]:
        get_vds_analysis_query: DocumentNode = gql(
            """
            query getVDSByAnalysisIds($vds_id: Int!) {
                analyses(id: {eq: $vds_id}) {
                    output
                    sequencingGroups {
                        id
                    }
                }
            }
        """,
        )
        query_results: dict[str, Any] = query(get_vds_analysis_query, variables={'vds_id': vds_id})
        vds_path: str = query_results['analyses'][0]['output']
        vds_sgids: list[str] = [sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']]
        return (vds_path, vds_sgids)

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        # Can't import it before all configs are set:
        from cpg_workflows.large_cohort import combiner

        workflow_config = config_retrieve('workflow')
        combiner_config = config_retrieve('combiner')

        output_vds_path: Path = self.expected_outputs(cohort)
        tmp_prefix = slugify(
            f"{self.tmp_prefix}/{workflow_config['cohort']}-{workflow_config['sequencing_type']}-{combiner_config['vds_version']}",
        )

        # create these as empty lists instead of None, they have the same truthiness
        vds_paths: list[str] = []
        sg_ids_in_vds: list[str] = []

        if combiner_config.get('vds_analysis_ids', None) is not None:
            for vds_id in combiner_config['vds_analysis_ids']:
                tmp_query_res, tmp_sg_ids_in_vds = self.get_vds_ids_output(vds_id)
                vds_paths.append(tmp_query_res)
                sg_ids_in_vds = sg_ids_in_vds + tmp_sg_ids_in_vds

        # Get SG IDs from the cohort object itself, rather than call Metamist.
        # Get VDS IDs first and filter out from this list
        cohort_sgs: list[SequencingGroup] = cohort.get_sequencing_groups(only_active=True)

        new_sg_gvcfs: list[str] = [str(sg.gvcf) for sg in cohort_sgs if sg.id not in sg_ids_in_vds]

        if len(new_sg_gvcfs) == 0 and len(vds_paths) <= 1:
            return self.make_outputs(cohort, self.expected_outputs(cohort))

        j: PythonJob = get_batch().new_python_job('Combiner', (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY})
        j.image(image_path('cpg_workflows'))
        j.memory(combiner_config['memory'])
        j.storage(combiner_config['storage'])

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(output_vds_path),
            sequencing_type=workflow_config['sequencing_type'],
            tmp_prefix=tmp_prefix,
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,  # this is a list or None, and in new_combiner None is made into []
            vds_paths=vds_paths,
        )

        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class SampleQC(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version() / 'sample_qc.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import sample_qc

        j = get_batch().new_job(
            'Sample QC',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                sample_qc,
                sample_qc.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(self.expected_outputs(cohort)),
                str(self.tmp_prefix),
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class DenseSubset(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return cohort.analysis_dataset.prefix() / get_workflow().name / dense_subset_version() / 'dense_subset.mt'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import dense_subset

        j = get_batch().new_job(
            'Dense Subset',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                dense_subset,
                dense_subset.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=DenseSubset, analysis_keys=['relatedness_ht'], analysis_type='custom')
class RelatednessPCRelate(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'relatedness_ht': (
                cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version() / 'relatedness.ht'
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(cohort)
        dense_mt_path = str(inputs.as_path(cohort, DenseSubset))
        dense_mt_name = dense_mt_path.split('/')[-1]

        # Estimate required storage - we create a new HT, so provision double the storage
        t_shirt_size_value = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        )
        required_storage = f'{t_shirt_size_value * 2}Gi'

        job = get_batch().new_job(f'Run Relatedness PCRelate stage: {cohort.name}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command('set -eux pipefail')

        # Localise the Dense MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {dense_mt_path} $BATCH_TMPDIR')

        required_cpu: int = config_retrieve(['RelatednessPCRelate', 'cores'], 8)
        job.cpu(required_cpu).storage(required_storage).memory('highmem')

        job.command(
            f'relatedness_pcrelate '
            f'--dense_mt "${{BATCH_TMPDIR}}/{dense_mt_name}" '
            f'--out {job.output} '
            f'--checkpoint "${{BATCH_TMPDIR}}"',
        )

        # Delocalise the output HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.output} {str(outputs["relatedness_ht"])}')

        return self.make_outputs(cohort, outputs, job)


@stage(
    required_stages=[SampleQC, RelatednessPCRelate],
    analysis_keys=['relateds_to_drop_ht'],
    analysis_type='custom',
)
class RelatednessFlag(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'relateds_to_drop_ht': (
                cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version() / 'relateds_to_drop.ht'
            ),
            'relateds_to_drop_gcta': (
                cohort.analysis_dataset.prefix()
                / get_workflow().name
                / relatedness_version()
                / 'gcta_relateds_remove.indi.list'
            ),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(cohort)

        t_shirt_size_value = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        )
        required_storage = f'{t_shirt_size_value * 2}Gi'

        job = get_batch().new_job(f'Run Relatedness Sample Flagging stage: {cohort.name}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command('set -eux pipefail')

        # Localise the Relatedness HT
        prcrelate_mt_path = str(inputs.as_path(cohort, RelatednessPCRelate, 'relatedness_ht'))
        prcrelate_mt_name = prcrelate_mt_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {prcrelate_mt_path} $BATCH_TMPDIR')

        sample_qc_ht_path = str(inputs.as_path(cohort, SampleQC))
        sample_qc_ht_name = sample_qc_ht_path.split('/')[-1]

        # Localise the SampleQC HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {sample_qc_ht_path} $BATCH_TMPDIR')

        required_cpu: int = config_retrieve(['RelatednessFlag', 'cores'], 8)
        job.cpu(required_cpu).storage(required_storage).memory('highmem')

        job.command(
            'relatedness_flag '
            f'--relatedness "${{BATCH_TMPDIR}}/{prcrelate_mt_name}" '
            f'--qc "${{BATCH_TMPDIR}}/{sample_qc_ht_name}" '
            f'--ht-out {job.relateds_to_drop_ht} '
            f'--gcta-out {job.relateds_to_drop_gcta} '
            f'--checkpoint "${{BATCH_TMPDIR}}" ',
        )

        # Delocalise the output HT
        job.command(
            f'gcloud --no-user-output-enabled storage cp -r {job.relateds_to_drop_ht} {str(outputs["relateds_to_drop_ht"])}',
        )
        job.command(
            f'gcloud --no-user-output-enabled storage cp -r {job.relateds_to_drop_gcta} {str(outputs["relateds_to_drop_gcta"])}',
        )

        return self.make_outputs(cohort, outputs, job)


@stage()
class DenseBackground(CohortStage):
    """
    Will densify a background dataset such as 1KG, HGDP, or 1KG-HGDP.
    This stage is run infrequently and only when the background dataset is updated.
    The densified background dataset is saved to that dataset's bucket as 'dense.mt'.
    Users can specify multiple background datasets to densify which will be saved to their respective buckets.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        background_datasets: list[str] = config_retrieve(['large_cohort', 'pca_background', 'datasets'], [])
        background_storage_paths: dict[str, Path] = {
            dataset: to_path(
                config_retrieve(['storage', dataset, 'default']),
            )
            for dataset in background_datasets
        }
        return {
            bg_dataset: background_storage_paths[bg_dataset] / 'dense' / 'dense.mt'
            for bg_dataset in background_datasets
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        if not (background_datasets := config_retrieve(['large_cohort', 'pca_background', 'datasets'], False)):
            return self.make_outputs(target=cohort)

        input_paths = {
            dataset: to_path(
                config_retrieve(['storage', dataset, 'default']),
            )
            / 'vds'
            # TODO: Allow for input background datasets to be '.mt' files
            / f'{dense_background_vds_version(dataset)}.vds'
            for dataset in background_datasets
        }

        job = get_batch().new_job(f'Densify {", ".join(background_datasets)} background datasets')
        job.storage('500Gi').image(config_retrieve(['workflow', 'driver_image']))

        outputs = self.expected_outputs(cohort)

        # Localise QC variants table
        qc_variants_ht = config_retrieve(['references', 'ancestry', 'sites_table'])

        jobs = []
        for bg_dataset in background_datasets:
            # Localise the background VDS
            job.command(
                f'gcloud --no-user-output-enabled storage cp -r {str(input_paths[bg_dataset])} $BATCH_TMPDIR',
            )
            job.command(
                'densify_background_dataset '
                f'--background-vds {str(input_paths[bg_dataset])} '
                f'--qc-variants-table {qc_variants_ht} '
                f'--dense-out {str(outputs[bg_dataset])} '
                f'--tmp {str(self.tmp_prefix)} ',
            )
            jobs.append(job)

        return self.make_outputs(target=cohort, jobs=jobs, data=outputs)


@stage(required_stages=[DenseSubset, SampleQC, DenseBackground])
class AncestryAddBackground(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'bg_qc_ht': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / sample_qc_version()
            / 'bg_sample_qc.ht',
            'bg_dense_mt': cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version() / 'bg_dense.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        dense_mt_path = inputs.as_path(cohort, DenseSubset)
        sample_qc_ht_path = inputs.as_path(cohort, SampleQC)
        dense_background_mt_path_dict: dict = inputs.as_dict(cohort, DenseBackground)

        outputs = self.expected_outputs(cohort)

        # Construct the background dataset parameters
        background_datasets = ' '.join(
            f'--background_dataset {name}={str(path)}' for name, path in dense_background_mt_path_dict.items()
        )

        job = get_batch().new_job('Add Ancestry Background')
        job.storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))
        job.command(
            f'ancestry_add_background '
            f'--qc_in {str(sample_qc_ht_path)} '
            f'--dense_in {str(dense_mt_path)} '
            f'--qc_out {str(outputs["bg_qc_ht"])} '
            f'--dense_out {str(outputs["bg_dense_mt"])} '
            f'{background_datasets} '
            f'--tmp {str(self.tmp_prefix)} ',
        )

        return self.make_outputs(target=cohort, jobs=job, data=outputs)


@stage(required_stages=[DenseSubset, RelatednessFlag])
class MakePlink(CohortStage):
    # Need to add AncestryBackground to the required stages only if we're using it.
    # Can't have this as default in decorator as it will trigger
    # the stage to run even if not required
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            self.required_stages_classes.append(AncestryAddBackground)

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return dict(
            root=get_workflow().prefix / 'gcta_pca' / 'PLINK' / gcta_version(),
            bed=get_workflow().prefix / 'gcta_pca' / 'PLINK' / f'{gcta_version()}.bed',
            bim=get_workflow().prefix / 'gcta_pca' / 'PLINK' / f'{gcta_version()}.bim',
            fam=get_workflow().prefix / 'gcta_pca' / 'PLINK' / f'{gcta_version()}.fam',
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        from cpg_workflows.large_cohort import make_plink

        j = get_batch().new_job('Make Plink Files', (self.get_job_attrs() or {}) | {'tool': 'hail query'})

        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                make_plink,
                make_plink.export_plink.__name__,
                str(inputs.as_path(cohort, DenseSubset)),
                str(self.expected_outputs(cohort)['root']),
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[MakePlink])
class GctaGRM(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return dict(
            grm_bin=get_workflow().prefix / 'gcta_pca' / 'GRM' / f'{gcta_version()}.grm.bin',
            grm_id=get_workflow().prefix / 'gcta_pca' / 'GRM' / f'{gcta_version()}.grm.id',
            grm_N_bin=get_workflow().prefix / 'gcta_pca' / 'GRM' / f'{gcta_version()}.grm.N.bin',
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        from cpg_workflows.jobs import gcta_GRM

        create_GRM_j = gcta_GRM.create_GRM(
            b=get_batch(),
            bed_file_path=str(inputs.as_path(cohort, MakePlink, key='bed')),
            bim_file_path=str(inputs.as_path(cohort, MakePlink, key='bim')),
            fam_file_path=str(inputs.as_path(cohort, MakePlink, key='fam')),
            output_path=str(self.expected_outputs(cohort)['grm_bin'].parent),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[create_GRM_j])


@stage(required_stages=[DenseSubset, RelatednessFlag, GctaGRM])
class GctaPCA(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return dict(
            pca_dir=get_workflow().prefix / 'gcta_pca' / 'PCA' / gcta_version(),
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.jobs import gcta_PCA

        ######### LOCALISATION #########
        # can't pass batch's

        outputs = self.expected_outputs(cohort)

        t_shirt_size_value = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        )
        required_storage = f'{t_shirt_size_value * 2}Gi'

        run_PCA_j: BashJob = get_batch().new_job(f'Run GCTA PCA: {cohort.name}')
        run_PCA_j.image(image_path('gcta'))
        run_PCA_j.command('set -eux pipefail')

        # Localise the GRM files
        grm_dir = str(inputs.as_path(cohort, GctaGRM, 'grm_bin').parent)
        gcta_relateds_to_drop_file = str(inputs.as_path(cohort, RelatednessFlag, 'relateds_to_drop_gcta'))
        relateds_name = gcta_relateds_to_drop_file.split('/')[-1]
        run_PCA_j.command(f'gcloud --no-user-output-enabled storage cp -r "{grm_dir}*" $BATCH_TMPDIR')
        run_PCA_j.command(f'gcloud --no-user-output-enabled storage cp -r {gcta_relateds_to_drop_file} $BATCH_TMPDIR')

        required_cpu: int = config_retrieve(['GctaPCA', 'cores'], 8)
        run_PCA_j.storage(required_storage).memory('highmem')
        run_PCA_j.command('ls -l $BATCH_TMPDIR')
        run_PCA_j.command('ls -l $BATCH_TMPDIR/GRM')
        run_PCA_j.command(
            f'gcta --grm $BATCH_TMPDIR/GRM/{gcta_version()} --remove $BATCH_TMPDIR/{relateds_name} --pca 10 --out {run_PCA_j.output}',
        )
        # run_PCA_j.command(
        #     'gcta_pca '
        #     '--grm-directory "${BATCH_TMPDIR}" '
        #     f'--output-path {run_PCA_j.output} '
        #     f'--version {gcta_version()} '
        #     f'--n-pcs {config_retrieve(["large_cohort", "n_pcs"])} '
        #     f'--relateds-to-drop "${{BATCH_TMPDIR}}/{relateds_name}" ',
        # )
        run_PCA_j.command(f'ls -l $BATCH_TMPDIR/{run_PCA_j.output}')
        # Delocalise the output
        run_PCA_j.command(f'gcloud --no-user-output-enabled storage cp -r {run_PCA_j.output} {str(outputs["pca_dir"])}')

        #################################
        ######### QUERY COMMAND #########

        # j = get_batch().new_job('Run GCTA PCA', (self.get_job_attrs() or {}) | {'tool': 'hail query'})
        # j.image(image_path('cpg_workflows'))
        # j.command(
        #     query_command(
        #         gcta_PCA,
        #         gcta_PCA.run_PCA.__name__,
        #         get_batch(),
        #         str(inputs.as_path(cohort, GctaGRM, 'grm_dir')),
        #         str(self.expected_outputs(cohort)['pca_dir']),
        #         gcta_version(),
        #         config_retrieve(['large_cohort', 'n_pcs']),
        #         str(inputs.as_path(cohort, RelatednessFlag, 'relateds_to_drop')),
        #         setup_gcp=True,
        #     ),
        # )

        #################################
        ######### FUNCTION CALL #########
        # b = get_batch()
        # # grm_directory = str(inputs.as_path(cohort, GctaGRM, 'grm_bin').parent)
        # bfile = b.read_input_group(
        #     **{
        #         'grm.bin': str(inputs.as_path(cohort, GctaGRM, 'grm_bin')),
        #         'grm.id': str(inputs.as_path(cohort, GctaGRM, 'grm_id')),
        #         'grm.N.bin': str(inputs.as_path(cohort, GctaGRM, 'grm_N_bin')),
        #     },
        # )
        # grm_directory = {
        #     'grm.bin': str(inputs.as_path(cohort, GctaGRM, 'grm_bin')),
        #     'grm.id': str(inputs.as_path(cohort, GctaGRM, 'grm_id')),
        #     'grm.N.bin': str(inputs.as_path(cohort, GctaGRM, 'grm_N_bin')),
        # }
        # logging.info(
        #     f'GRM files LARGE COHORT: bin: {bfile["grm.bin"]}, id: {bfile["grm.id"]}, N: {bfile["grm.N.bin"]}',
        # )
        # run_PCA_j = gcta_PCA.run_PCA(
        #     b=b,
        #     grm_directory=grm_directory,
        #     output_path=str(self.expected_outputs(cohort)['pca_dir']),
        #     version=gcta_version(),
        #     n_pcs=config_retrieve(['large_cohort', 'n_pcs']),
        #     relateds_to_drop=str(inputs.as_path(cohort, RelatednessFlag, 'relateds_to_drop_gcta')),
        # )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[run_PCA_j])


@stage(required_stages=[SampleQC, RelatednessFlag, DenseSubset])
class AncestryPCA(CohortStage):
    # Need to add AncestryBackground to the required stages only if we're using it.
    # Can't have this as default in decorator as it will trigger
    # the stage to run even if not required
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            self.required_stages_classes.append(AncestryAddBackground)

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'scores': ancestry_prefix / 'scores.ht',
            'eigenvalues': ancestry_prefix / 'eigenvalues.ht',
            'loadings': ancestry_prefix / 'loadings.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        # Decide which stage to pull input from
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            dense_mt_path = str(inputs.as_path(cohort, AncestryAddBackground, 'bg_dense_mt'))
        else:
            dense_mt_path = str(inputs.as_path(cohort, DenseSubset))

        related = str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop_ht'))

        job = get_batch().new_bash_job('Run Ancestry PCA')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('50Gi').memory('highmem').cpu(config_retrieve(['Ancestry', 'cores'], 8))

        # Localise the dense MT
        dense_name = dense_mt_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {dense_mt_path} $BATCH_TMPDIR')

        # Localise the relatedness HT
        related_name = related.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {related} $BATCH_TMPDIR')

        job.command(
            'ancestry_pca '
            f'--dense_mt "${{BATCH_TMPDIR}}/{dense_name}" '
            f'--related "${{BATCH_TMPDIR}}/{related_name}" '
            f'--scores_out {job.scores} '
            f'--eigen_out {job.eigen} '
            f'--loadings_out {job.loadings} ',
        )

        # Copy 3 tables out
        outputs = self.expected_outputs(cohort)
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.scores} {str(outputs["scores"])}')
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.eigen} {str(outputs["eigenvalues"])}')
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.loadings} {str(outputs["loadings"])}')

        return self.make_outputs(target=cohort, data=outputs, jobs=job)


@stage(required_stages=[AncestryPCA, SampleQC])
class AncestryInfer(CohortStage):
    # Need to add AncestryBackground to the required stages only if we're using it.
    # Can't have this as default in decorator as it will trigger
    # the stage to run even if not required
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            self.required_stages_classes.append(AncestryAddBackground)

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:

        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'inferred_pop': ancestry_prefix / 'inferred_pop.ht',
            'sample_qc_ht': ancestry_prefix / 'sample_qc_ht.ht',
            'pickled_model': ancestry_prefix / 'pop.RF_fit.pickle',
            'inference_txt': ancestry_prefix / 'RF_pop_assignments.txt.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        if config_retrieve(['large_cohort', 'pca_background', 'datasets'], False):
            sample_qc_ht_path = str(inputs.as_path(cohort, AncestryAddBackground, 'bg_qc_ht'))
        else:
            sample_qc_ht_path = str(inputs.as_path(cohort, SampleQC))
        scores_table = str(inputs.as_path(cohort, AncestryPCA, 'scores'))

        job = get_batch().new_bash_job('Ancestry Infer Populations')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('50Gi').memory('highmem').cpu(config_retrieve(['Ancestry', 'cores'], 8))

        sample_qc_name = sample_qc_ht_path.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {sample_qc_ht_path} $BATCH_TMPDIR')

        scores_name = scores_table.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {scores_table} $BATCH_TMPDIR')

        job.command(
            'ancestry_infer_labels '
            f'--scores_ht "${{BATCH_TMPDIR}}/{scores_name}" '
            f'--qc_in "${{BATCH_TMPDIR}}/{sample_qc_name}" '
            f'--qc_out {job.qc} '
            f'--pop_ht_out {job.ht} '
            f'--pickle_out {job.pickle} '
            f'--txt_out {job.txt} ',
        )

        outputs = self.expected_outputs(cohort)

        # Delocalise the scores HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.ht} {str(outputs["inferred_pop"])}')

        # Delocalise the SampleQC HT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {job.qc} {str(outputs["sample_qc_ht"])}')

        get_batch().write_output(job.pickle, str(outputs['pickled_model']))
        get_batch().write_output(job.txt, str(outputs['inference_txt']))

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=[SampleQC, DenseSubset, RelatednessFlag])
class Ancestry(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        ancestry_prefix = get_workflow().prefix / 'ancestry'
        return {
            'scores': ancestry_prefix / 'scores.ht',
            'eigenvalues': ancestry_prefix / 'eigenvalues.ht',
            'loadings': ancestry_prefix / 'loadings.ht',
            'inferred_pop': ancestry_prefix / 'inferred_pop.ht',
            'sample_qc_ht': ancestry_prefix / 'sample_qc_ht.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_pca import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop_ht'),
                tmp_prefix=self.tmp_prefix,
                out_scores_ht_path=self.expected_outputs(cohort)['scores'],
                out_eigenvalues_ht_path=self.expected_outputs(cohort)['eigenvalues'],
                out_loadings_ht_path=self.expected_outputs(cohort)['loadings'],
                out_inferred_pop_ht_path=self.expected_outputs(cohort)['inferred_pop'],
                out_sample_qc_ht_path=self.expected_outputs(cohort)['sample_qc_ht'],
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, AncestryInfer, RelatednessFlag])
class AncestryPlots(CohortStage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.out_prefix = get_workflow().web_prefix / 'ancestry'
        self.out_fname_pattern = '{scope}_pc{pci}_{pca_suffix}.{ext}'

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        n_pcs = get_config()['large_cohort']['n_pcs']
        # if there is a pca_plot_name given, add this to the output name
        plot_name = get_config()['large_cohort'].get('pca_plot_name')
        pca_suffix = ''
        if plot_name:
            pca_suffix = plot_name.replace('-', '_')
        return {
            str(pc_num): self.out_prefix
            / self.out_fname_pattern.format(
                scope='dataset',
                pci=pc_num,
                pca_suffix=pca_suffix,
                ext='html',
            )
            for pc_num in range(1, n_pcs)
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_plots import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                out_path_pattern=self.out_prefix / self.out_fname_pattern,
                sample_qc_ht_path=inputs.as_path(cohort, Ancestry, key='sample_qc_ht'),
                scores_ht_path=inputs.as_path(cohort, Ancestry, key='scores'),
                eigenvalues_ht_path=inputs.as_path(cohort, Ancestry, key='eigenvalues'),
                loadings_ht_path=inputs.as_path(cohort, Ancestry, key='loadings'),
                inferred_pop_ht_path=inputs.as_path(
                    cohort,
                    Ancestry,
                    key='inferred_pop',
                ),
                relateds_to_drop_ht_path=inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop_ht'),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner, SampleQC, RelatednessFlag])
class MakeSiteOnlyVcf(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'siteonly.vcf.bgz',
            'tbi': self.tmp_prefix / 'siteonly.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import site_only_vcf

        j = get_batch().new_job(
            'MakeSiteOnlyVcf',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                site_only_vcf,
                site_only_vcf.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop_ht')),
                str(self.expected_outputs(cohort)['vcf']),
                str(self.tmp_prefix),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=MakeSiteOnlyVcf)
class Vqsr(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'vcf': self.tmp_prefix / 'siteonly.vqsr.vcf.gz',
            'tbi': self.tmp_prefix / 'siteonly.vqsr.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.jobs import vqsr

        vcf_path = inputs.as_path(cohort, MakeSiteOnlyVcf, key='vcf')
        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=vcf_path,
            gvcf_count=len(cohort.get_sequencing_groups()),
            out_path=self.expected_outputs(cohort)['vcf'],
            tmp_prefix=self.tmp_prefix,
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=Vqsr)
class LoadVqsr(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return get_workflow().prefix / 'vqsr.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import load_vqsr

        j = get_batch().new_job(
            'LoadVqsr',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                load_vqsr,
                load_vqsr.run.__name__,
                str(inputs.as_path(cohort, Vqsr, key='vcf')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[Combiner, SampleQC, RelatednessFlag])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return get_workflow().prefix / 'frequencies.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import frequencies

        j = get_batch().new_job(
            'Frequencies',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                frequencies,
                frequencies.run.__name__,
                str(inputs.as_path(cohort, Combiner)),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, RelatednessFlag, key='relateds_to_drop_ht')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])
