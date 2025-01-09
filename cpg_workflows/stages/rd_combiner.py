"""
All Stages relating to the seqr_loader pipeline, reimplemented from scratch to
use the gVCF combiner instead of joint-calling.
"""

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.gcloud_composer import gcloud_compose_vcf_from_manifest
from cpg_workflows.jobs.rd_combiner import combiner
from cpg_workflows.jobs.rd_combiner.vep import add_vep_jobs
from cpg_workflows.jobs.rd_combiner.vqsr import (
    apply_recalibration_indels,
    apply_snp_vqsr_to_fragments,
    gather_tranches,
    train_vqsr_indels,
    train_vqsr_snp_tranches,
    train_vqsr_snps,
)
from cpg_workflows.targets import Dataset, MultiCohort
from cpg_workflows.utils import get_all_fragments_from_manifest, get_logger
from cpg_workflows.workflow import (
    DatasetStage,
    MultiCohortStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

LATEST_ANALYSIS_QUERY = gql(
    """
    query LatestAnalysisEntry($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                meta
                output
                sequencingGroups {
                    id
                }
                timestampCompleted
            }
        }
    }
""",
)
SPECIFIC_VDS_QUERY = gql(
    """
    query getVDSByAnalysisId($vds_id: Int!) {
        analyses(id: {eq: $vds_id}) {
            output
            sequencingGroups {
                id
            }
        }
    }
""",
)
SHARD_MANIFEST = 'shard-manifest.txt'


def query_for_specific_vds(vds_id: int) -> tuple[str, set[str]] | None:
    """
    query for a specific analysis of type entry_type for a dataset
    if found, return the set of SG IDs in the VDS (using the metadata)

    - stolen from the cpg_workflows.large_cohort.combiner Stage, but duplicated here so we can split pipelines without
      further code changes

    Args:
        vds_id (int): analysis id to query for

    Returns:
        either None if the analysis wasn't found, or a set of SG IDs in the VDS
    """

    # query for the exact, single analysis entry
    query_results: dict[str, dict] = query(SPECIFIC_VDS_QUERY, variables={'vds_id': vds_id})

    if not query_results['analyses']:
        return None
    vds_path: str = query_results['analyses'][0]['output']
    sg_ids = {sg['id'] for sg in query_results['analyses'][0]['sequencingGroups']}
    return vds_path, sg_ids


def query_for_latest_vds(dataset: str, entry_type: str = 'combiner') -> dict | None:
    """
    query for the latest analysis of type entry_type for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest analysis
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(LATEST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': entry_type})

    analyses_by_date = {}

    for analysis in result['project']['analyses']:
        if analysis['output'] and (
            analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type'])
        ):
            analyses_by_date[analysis['timestampCompleted']] = analysis

    if not analyses_by_date:
        get_logger(__file__).warning(f'No analysis of type {entry_type} found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort as strings
    return analyses_by_date[sorted(analyses_by_date)[-1]]


@stage(analysis_type='combiner', analysis_keys=['vds'])
class CreateVdsFromGvcfsWithHailCombiner(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'combiner_plan': str(self.tmp_prefix / 'combiner_plan.json'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        outputs: dict[str, str | Path] = self.expected_outputs(multicohort)

        # create these as empty lists instead of None, they have the same truthiness
        vds_path: str | None = None
        sg_ids_in_vds: set[str] = set()

        # check for a VDS by ID
        if vds_id := config_retrieve(['workflow', 'use_specific_vds'], False):
            vds_result_or_none = query_for_specific_vds(vds_id)
            if vds_result_or_none is None:
                raise ValueError(f'Specified VDS ID {vds_id} not found in Metamist')

            # if not none, unpack the result
            vds_path, sg_ids_in_vds = vds_result_or_none

        # check for existing VDS by getting all and fetching latest
        elif config_retrieve(['workflow', 'check_for_existing_vds'], True):
            get_logger(__file__).info('Checking for existing VDS')
            if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
                vds_path = existing_vds_analysis_entry['output']
                sg_ids_in_vds = {sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']}

        else:
            get_logger(__file__).info('Not continuing from any previous VDS, creating new Combiner from gVCFs only')

        new_sg_gvcfs: list[str] = [
            str(sg.gvcf)
            for sg in multicohort.get_sequencing_groups()
            if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        ]

        if not new_sg_gvcfs:
            get_logger(__file__).info('No GVCFs to combine')
            get_logger(__file__).info(f'Checking if VDS exists: {outputs["vds"]}: {outputs["vds"].exists()}')  # type: ignore
            return self.make_outputs(multicohort, outputs)

        combiner_job = get_batch().new_python_job('CreateVdsFromGvcfsWithHailCombiner', {'stage': self.name})
        combiner_job.image(config_retrieve(['workflow', 'driver_image']))
        combiner_job.memory(config_retrieve(['combiner', 'driver_memory'], 'highmem'))
        combiner_job.storage(config_retrieve(['combiner', 'driver_storage']))
        combiner_job.cpu(config_retrieve(['combiner', 'driver_cores'], 2))

        # Default to GRCh38 for reference if not specified
        combiner_job.call(
            combiner.run,
            output_vds_path=str(outputs['vds']),
            save_path=outputs['combiner_plan'],
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_paths=[vds_path] if vds_path else None,
            force_new_combiner=config_retrieve(['combiner', 'force_new_combiner'], False),
        )

        return self.make_outputs(multicohort, outputs, combiner_job)


@stage(required_stages=[CreateVdsFromGvcfsWithHailCombiner], analysis_type='matrixtable', analysis_keys=['mt'])
class CreateDenseMtFromVdsWithHail(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        the MT and both shard_manifest files are Paths, so this stage will rerun if any of those are missing
        the VCFs are written as a directory, rather than a single VCF, so we can't check its existence well

        Needs a range of INFO fields to be present in the VCF
        """
        prefix = self.prefix
        temp_prefix = self.tmp_prefix

        return {
            'mt': prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF (header-per-shard)
            'hps_vcf_dir': str(prefix / f'{multicohort.name}.vcf.bgz'),
            # this will be the file which contains the name of all fragments (header-per-shard)
            'hps_shard_manifest': prefix / f'{multicohort.name}.vcf.bgz' / SHARD_MANIFEST,
            # this will be the write path for fragments of sites-only VCF (separate header)
            'separate_header_vcf_dir': str(temp_prefix / f'{multicohort.name}_separate.vcf.bgz'),
            # this will be the file which contains the name of all fragments (separate header)
            'separate_header_manifest': temp_prefix / f'{multicohort.name}_separate.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        output = self.expected_outputs(multicohort)

        # partitions to coalesce the data into
        partitions = config_retrieve(['workflow', 'densify_partitions'], 2500)

        # not currently in use (see #1078)
        partition_strategy = config_retrieve(['workflow', 'partition_strategy'], 'naive')

        densify_job = get_batch().new_job('CreateDenseMtFromVdsWithHail')
        densify_job.image(config_retrieve(['workflow', 'driver_image']))
        densify_job.command(
            'mt_from_vds '
            f'--input {str(inputs.as_dict(multicohort, CreateVdsFromGvcfsWithHailCombiner)["vds"])} '
            f'--output {str(output["mt"])} '
            f'--partitions {partitions} '
            f'--partition_strategy {partition_strategy} '
            f'--sites_only {output["hps_vcf_dir"]} '
            f'--separate_header {output["separate_header_vcf_dir"]} ',
        )
        return self.make_outputs(multicohort, output, densify_job)


@stage(analysis_keys=['vcf'], required_stages=[CreateDenseMtFromVdsWithHail])
class ConcatenateVcfFragmentsWithGcloud(MultiCohortStage):
    """
    Takes a manifest of VCF fragments, and produces a single VCF file
    This is disconnected from the previous stage, but requires it to be run first
    So we check for the exact same output, and fail if we're not ready to start
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.tmp_prefix / 'gcloud_composed_sitesonly.vcf.bgz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs to take a manifest of VCF fragments, and produce a single VCF file
        The VCF being composed here has a single header in a separate file, the first entry in the manifest
        This means we can combine the VCF header and data fragments through concatenation
        and the result will be a spec-compliant VCF
        """
        manifest_file = inputs.as_path(
            target=multicohort,
            stage=CreateDenseMtFromVdsWithHail,
            key='separate_header_manifest',
        )
        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the rd_combiner workflow with workflows.last_stages=[CreateDenseMtFromVdsWithHail]',
            )

        outputs = self.expected_outputs(multicohort)

        jobs = gcloud_compose_vcf_from_manifest(
            manifest_path=manifest_file,
            intermediates_path=str(self.tmp_prefix / 'temporary_compose_intermediates'),
            output_path=str(outputs),
            job_attrs={'stage': self.name},
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(required_stages=[ConcatenateVcfFragmentsWithGcloud])
class TrainVqsrIndelModelOnCombinerData(MultiCohortStage):
    """
    Train VQSR Indel model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path | str]:
        prefix = self.prefix
        return {
            'indel_recalibrations': prefix / 'indel.recal',
            'indel_tranches': prefix / 'indel.tranches',
            'indel_prefix': str(prefix / 'indel'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_path(multicohort, ConcatenateVcfFragmentsWithGcloud)
        outputs = self.expected_outputs(multicohort)
        indel_calibration_job = train_vqsr_indels(
            sites_only_vcf=str(composed_sitesonly_vcf),
            output_prefix=str(outputs['indel_prefix']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=indel_calibration_job)


@stage(required_stages=[ConcatenateVcfFragmentsWithGcloud])
class TrainVqsrSnpModelOnCombinerData(MultiCohortStage):
    """
    Train VQSR SNP model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.prefix / 'snp_model'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_path(multicohort, ConcatenateVcfFragmentsWithGcloud)
        outputs = self.expected_outputs(multicohort)
        snp_calibration_job = train_vqsr_snps(
            sites_only_vcf=str(composed_sitesonly_vcf),
            snp_model=str(outputs),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=snp_calibration_job)


@stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpModelOnCombinerData])
class TrainVqsrSnpTranches(MultiCohortStage):
    """
    Scattered training of VQSR tranches for SNPs
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, str | Path]:

        prefix = self.tmp_prefix
        return {
            'tranche_marker': prefix / 'tranches_trained',
            'temp_path': str(prefix / 'vqsr_snp_tranches'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')
        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the rd_combiner workflow with workflows.last_stages=[CreateDenseMtFromVdsWithHail]',
            )

        outputs = self.expected_outputs(multicohort)
        snp_model_path = inputs.as_path(target=multicohort, stage=TrainVqsrSnpModelOnCombinerData)

        jobs = train_vqsr_snp_tranches(
            manifest_file=manifest_file,
            snp_model_path=str(snp_model_path),
            output_path=str(outputs['tranche_marker']),
            temp_path=to_path(outputs['temp_path']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(required_stages=[CreateDenseMtFromVdsWithHail, TrainVqsrSnpTranches])
class GatherTrainedVqsrSnpTranches(MultiCohortStage):
    """
    Scattered training of VQSR tranches for SNPs
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.prefix / 'snp_tranches'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')
        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the rd_combiner workflow with workflows.last_stages=[CreateDenseMtFromVdsWithHail]',
            )

        outputs = self.expected_outputs(multicohort)

        jobs = gather_tranches(
            manifest_file=manifest_file,
            temp_path=to_path(inputs.as_str(target=multicohort, stage=TrainVqsrSnpTranches, key='temp_path')),
            output_path=str(outputs),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        GatherTrainedVqsrSnpTranches,
        TrainVqsrSnpTranches,
    ],
)
class RunTrainedSnpVqsrOnCombinerFragments(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.tmp_prefix / 'vqsr.vcf.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')
        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the rd_combiner workflow with workflows.last_stages=[CreateDenseMtFromVdsWithHail]',
            )

        outputs = self.expected_outputs(multicohort)
        tranche_recal_temp = to_path(inputs.as_str(target=multicohort, stage=TrainVqsrSnpTranches, key='temp_path'))
        tranche_file = inputs.as_path(target=multicohort, stage=GatherTrainedVqsrSnpTranches)

        jobs = apply_snp_vqsr_to_fragments(
            manifest_file=manifest_file,
            tranche_file=str(tranche_file),
            temp_path=tranche_recal_temp,
            output_path=str(outputs),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(
    analysis_type='qc',
    required_stages=[
        RunTrainedSnpVqsrOnCombinerFragments,
        TrainVqsrIndelModelOnCombinerData,
    ],
)
class RunTrainedIndelVqsrOnCombinedVcf(MultiCohortStage):
    """
    Run Indel VQSR on the reconstituted, SNP-annotated, VCF
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return self.prefix / 'vqsr_snps_and_indels.vcf.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        annotated_vcf = inputs.as_path(target=multicohort, stage=RunTrainedSnpVqsrOnCombinerFragments)
        indel_recalibrations = inputs.as_path(
            target=multicohort,
            stage=TrainVqsrIndelModelOnCombinerData,
            key='indel_recalibrations',
        )
        indel_tranches = inputs.as_path(
            target=multicohort,
            stage=TrainVqsrIndelModelOnCombinerData,
            key='indel_tranches',
        )

        indel_recal_job = apply_recalibration_indels(
            snp_annotated_vcf=annotated_vcf,
            indel_recalibration=indel_recalibrations,
            indel_tranches=indel_tranches,
            output_path=outputs,
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=indel_recal_job)


@stage(analysis_type='custom', required_stages=[CreateDenseMtFromVdsWithHail])
class AnnotateFragmentedVcfWithVep(MultiCohortStage):
    """
    Annotate VCF with VEP.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        """
        Should this be in tmp? We'll never use it again maybe?
        """
        return self.prefix / 'vep.ht'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        manifest_file = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='hps_shard_manifest')
        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the rd_combiner workflow with workflows.last_stages=[CreateDenseMtFromVdsWithHail]',
            )

        input_vcfs = get_all_fragments_from_manifest(manifest_file)

        if len(input_vcfs) == 0:
            raise ValueError(f'No VCFs in {manifest_file}')

        vep_jobs = add_vep_jobs(
            input_vcfs=input_vcfs,
            final_out_path=outputs,
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=vep_jobs)


@stage(
    analysis_type='matrixtable',
    required_stages=[
        CreateDenseMtFromVdsWithHail,
        AnnotateFragmentedVcfWithVep,
        RunTrainedIndelVqsrOnCombinedVcf,
    ],
)
class AnnotateCohortSmallVariantsWithHailQuery(MultiCohortStage):
    """
    Annotate small variants with VEP and VQSR
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        """
        Expected to write a matrix table.
        """
        return self.tmp_prefix / 'annotate_cohort.mt'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """

        Args:
            multicohort ():
            inputs ():
        """

        outputs = self.expected_outputs(multicohort)
        vep_ht_path = inputs.as_path(target=multicohort, stage=AnnotateFragmentedVcfWithVep)
        vqsr_vcf = inputs.as_path(target=multicohort, stage=RunTrainedIndelVqsrOnCombinedVcf)
        variant_mt = inputs.as_path(target=multicohort, stage=CreateDenseMtFromVdsWithHail, key='mt')

        job = get_batch().new_job(self.name, self.get_job_attrs(multicohort))
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2).memory('highmem').storage('10Gi')
        job.command(
            'annotate_cohort_small '
            f'--input {variant_mt} '
            f'--output {outputs} '
            f'--vep {vep_ht_path} '
            f'--checkpoint {self.tmp_prefix} '
            f'--vqsr {vqsr_vcf} ',
        )
        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[AnnotateCohortSmallVariantsWithHailQuery], analysis_type='matrixtable', analysis_keys=['mt'])
class SubsetMatrixTableToDatasetUsingHailQuery(DatasetStage):
    """
    Subset the MT to a single dataset
    Skips this stage if the MultiCohort has only one dataset
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path] | None:
        """
        Expected to generate a MatrixTable
        This is kinda transient, so shove it in tmp
        """
        if len(get_multicohort().get_datasets()) == 1:
            get_logger().info(f'Skipping SubsetMatrixTableToDataset for single Dataset {dataset}')
            return None
        return {
            'mt': dataset.tmp_prefix() / 'mt' / self.name / f'{get_workflow().output_version}-{dataset.name}.mt',
            'id_file': dataset.tmp_prefix() / f'{get_workflow().output_version}-{dataset.name}-SG-ids.txt',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        outputs = self.expected_outputs(dataset)

        # only create dataset MTs for datasets specified in the config
        # and only run this stage if the callset has multiple datasets
        if (outputs is None) or (
            dataset.name not in config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
        ):
            get_logger().info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        variant_mt = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSmallVariantsWithHailQuery)

        # write a list of all the SG IDs to retain
        if not config_retrieve(['workflow', 'dry_run'], False):
            with outputs['id_file'].open('w') as f:
                for sg in dataset.get_sequencing_groups():
                    f.write(f'{sg.id}\n')

        job = get_batch().new_job(self.name, self.get_job_attrs(dataset))
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2).memory('highmem').storage('10Gi')
        job.command(
            'subset_mt_to_dataset '
            f'--input {str(variant_mt)} '
            f'--output {str(outputs["mt"])} '
            f'--sg_id_file {str(outputs["id_file"])} ',
        )
        return self.make_outputs(dataset, data=outputs, jobs=job)


@stage(
    required_stages=[
        AnnotateCohortSmallVariantsWithHailQuery,
        SubsetMatrixTableToDatasetUsingHailQuery,
    ],
    analysis_type='matrixtable',
)
class AnnotateDatasetSmallVariantsWithHailQuery(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        return dataset.prefix() / 'mt' / self.name / f'{get_workflow().output_version}-{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        # only create final MTs for datasets specified in the config
        if dataset.name not in config_retrieve(['workflow', 'write_mt_for_datasets'], default=[]):
            get_logger().info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        # choose the input MT based on the number of datasets in the MultiCohort
        if len(get_multicohort().get_datasets()) == 1:
            input_mt = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSmallVariantsWithHailQuery)
        else:
            input_mt = inputs.as_path(target=dataset, stage=SubsetMatrixTableToDatasetUsingHailQuery, key='mt')

        outputs = self.expected_outputs(dataset)

        job = get_batch().new_job(self.name, self.get_job_attrs(dataset))
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2).memory('highmem').storage('10Gi')
        job.command(f'annotate_dataset_small --input {str(input_mt)} --output {str(outputs)} ')
        return self.make_outputs(dataset, data=outputs, jobs=job)
