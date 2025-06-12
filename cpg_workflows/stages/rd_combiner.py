"""
All Stages relating to the seqr_loader pipeline, reimplemented from scratch to
use the gVCF combiner instead of joint-calling.
"""

from functools import cache

from google.api_core.exceptions import PermissionDenied

from hail.vds import read_vds

from cpg_utils import Path, to_path
from cpg_utils.cloud import read_secret
from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import get_batch, init_batch
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
from cpg_workflows.jobs.seqr_loader import cohort_to_vcf_job
from cpg_workflows.targets import Dataset, MultiCohort
from cpg_workflows.utils import exists, get_all_fragments_from_manifest, get_logger, tshirt_mt_sizing
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


def manually_find_ids_from_vds(vds_path: str) -> set[str]:
    """
    during development and the transition to input_cohorts over input_datasets, there are some instances
    where we have VDS entries in Metamist, but the analysis entry contains SG IDs which weren't combined into the VDS

    This check bypasses the quick "get all SG IDs in the VDS analysis entry" check,
    and instead checks the exact contents of the VDS

    Args:
        vds_path (str): path to the VDS. Assuming it exists, this will be checked before calling this method

    Returns:
        set[str]: the set of sample IDs in the VDS
    """
    init_batch()
    vds = read_vds(vds_path)

    # find the samples in the Variant Data MT
    return set(vds.variant_data.s.collect())


@cache
def get_family_sequencing_groups(dataset: Dataset) -> dict | None:
    """
    Get the subset of sequencing groups that are in the specified families for a dataset
    Returns a dict containing the sequencing groups and a name suffix for the outputs
    """
    if not config_retrieve(['workflow', dataset.name, 'only_families'], []):
        return None
    only_family_ids = set(config_retrieve(['workflow', dataset.name, 'only_families'], []))
    # keep only the SG IDs for the families in the only_families list
    get_logger().info(f'Finding sequencing groups for families {only_family_ids} in dataset {dataset.name}')
    family_sg_ids = [sg.id for sg in dataset.get_sequencing_groups() if sg.pedigree.fam_id in only_family_ids]
    if not family_sg_ids:
        raise ValueError(f'No sequencing groups found for families {only_family_ids} in dataset {dataset.name}.')
    get_logger().info(f'Keeping only {len(family_sg_ids)} SGs from families {len(only_family_ids)} in {dataset}:')
    get_logger().info(only_family_ids)
    get_logger().info(family_sg_ids)

    import hashlib

    h = hashlib.sha256(''.join(sorted(family_sg_ids)).encode()).hexdigest()[:4]
    name_suffix = f'{len(family_sg_ids)}_sgs-{len(only_family_ids)}_families-{h}'

    return {'family_sg_ids': family_sg_ids, 'name_suffix': name_suffix}


@stage(analysis_type='combiner', analysis_keys=['vds'])
class CreateVdsFromGvcfsWithHailCombiner(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'combiner_plan': str(self.tmp_prefix / 'combiner_plan.json'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        outputs = self.expected_outputs(multicohort)

        # we only ever build on top of a single VDS, or start from scratch
        vds_path: str | None = None

        # create these as empty iterables
        sg_ids_in_vds: set[str] = set()
        sgs_to_remove: list[str] = []

        # use a VDS path from config file, if possible
        if vds_path := config_retrieve(['workflow', 'specific_vds'], None):
            get_logger().info(f'Using VDS path from config: {vds_path}')

        # check for existing VDS by getting all and fetching latest
        elif config_retrieve(['workflow', 'check_for_existing_vds']):
            get_logger().info('Checking for existing VDS')
            if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
                vds_path = existing_vds_analysis_entry['output']
                sg_ids_in_vds = {sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']}

        else:
            get_logger().info('Not continuing from any previous VDS, creating new Combiner from gVCFs only')

        # quick check - if we found a VDS, guarantee it exists
        if vds_path and not exists(vds_path):
            raise ValueError(f'VDS {vds_path} does not exist, but has an Analysis Entry')

        # this is a quick and confident check on current VDS contents, but does require a direct connection to the VDS
        # by default this is True, and can be deactivated in config
        if vds_path and config_retrieve(['workflow', 'manually_check_vds_sg_ids']):
            sg_ids_in_vds = manually_find_ids_from_vds(vds_path)

        # technicality; this can be empty - in a situation where we have a VDS already and the current MCohort has FEWER
        # SG IDs in it, this stage will re-run because the specific hash will be different to the previous VDS
        # See https://github.com/populationgenomics/production-pipelines/issues/1126
        new_sg_gvcfs: list[str] = [
            str(sg.gvcf)
            for sg in multicohort.get_sequencing_groups()
            if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        ]

        # final check - if we have a VDS, and we have a current MultiCohort
        # detect any samples which should be _removed_ from the current VDS prior to further combining taking place
        if sg_ids_in_vds:
            sgs_in_mc: list[str] = multicohort.get_sequencing_group_ids()
            get_logger().info(f'Found {len(sg_ids_in_vds)} SG IDs in VDS {vds_path}')
            get_logger().info(f'Total {len(sgs_in_mc)} SGs in this MultiCohort')

            sgs_to_remove = sorted(set(sg_ids_in_vds) - set(sgs_in_mc))

            if sgs_to_remove:
                get_logger().info(f'Removing {len(sgs_to_remove)} SGs from VDS {vds_path}')
                get_logger().info(f'SGs to remove: {sgs_to_remove}')

        if not (new_sg_gvcfs or sgs_to_remove):
            get_logger().info('No GVCFs to add to, or remove from, existing VDS')
            get_logger().info(f'Checking if VDS exists: {outputs["vds"]}: {outputs["vds"].exists()}')  # type: ignore
            return self.make_outputs(multicohort, outputs)

        combiner_job = get_batch().new_python_job('CreateVdsFromGvcfsWithHailCombiner', {'stage': self.name})
        combiner_job.image(config_retrieve(['workflow', 'driver_image']))
        combiner_job.memory(config_retrieve(['combiner', 'driver_memory']))
        combiner_job.storage(config_retrieve(['combiner', 'driver_storage']))
        combiner_job.cpu(config_retrieve(['combiner', 'driver_cores']))

        # set this job to be non-spot (i.e. non-preemptible)
        # previous issues with preemptible VMs led to multiple simultaneous QOB groups processing the same data
        combiner_job.spot(config_retrieve(['combiner', 'preemptible_vms']))

        # Default to GRCh38 for reference if not specified
        combiner_job.call(
            combiner.run,
            output_vds_path=str(outputs['vds']),
            save_path=outputs['combiner_plan'],
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_path=vds_path,
            force_new_combiner=config_retrieve(['combiner', 'force_new_combiner'], False),
            sgs_to_remove=sgs_to_remove,
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


@stage(required_stages=[AnnotateCohortSmallVariantsWithHailQuery])
class SubsetMatrixTableToDatasetUsingHailQuery(DatasetStage):
    """
    Subset the MT to a single dataset - or a subset of families within a dataset
    Skips this stage if the MultiCohort has only one dataset
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path] | None:
        """
        Expected to generate a MatrixTable
        This is kinda transient, so shove it in tmp.

        If subsetting to certain families, the output will be named accordingly
        """
        if family_sgs := get_family_sequencing_groups(dataset):
            return {
                'mt': dataset.tmp_prefix()
                / 'mt'
                / self.name
                / f'{get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.mt',
                'id_file': dataset.tmp_prefix()
                / f'{get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}-SG-ids.txt',
            }
        elif len(get_multicohort().get_datasets()) == 1:
            get_logger().info(f'Skipping SubsetMatrixTableToDataset for single Dataset {dataset}')
            return None
        else:
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

        if family_sgs := get_family_sequencing_groups(dataset):
            family_sg_ids = family_sgs['family_sg_ids']
        else:
            family_sg_ids = None

        # write a list of all the SG IDs to retain
        # don't try and extract samples which didn't have a gVCF
        if not config_retrieve(['workflow', 'dry_run'], False):
            with outputs['id_file'].open('w') as f:
                for sg in dataset.get_sequencing_groups():
                    if ((family_sg_ids is None) or sg.id in family_sg_ids) and (sg.gvcf is not None):
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
        if family_sgs := get_family_sequencing_groups(dataset):
            return (
                dataset.prefix()
                / 'mt'
                / self.name
                / f'{get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.mt'
            )
        else:
            return dataset.prefix() / 'mt' / self.name / f'{get_workflow().output_version}-{dataset.name}.mt'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        # only create final MTs for datasets specified in the config
        if dataset.name not in config_retrieve(['workflow', 'write_mt_for_datasets'], default=[]):
            get_logger().info(f'Skipping AnnotateDataset mt subsetting for {dataset}')
            return self.make_outputs(dataset)

        family_sgs = get_family_sequencing_groups(dataset)
        # choose the input MT based on the number of datasets in the MultiCohort and the presence of family SGs
        if len(get_multicohort().get_datasets()) == 1 and family_sgs is None:
            input_mt = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSmallVariantsWithHailQuery)
        else:
            input_mt = inputs.as_path(target=dataset, stage=SubsetMatrixTableToDatasetUsingHailQuery, key='mt')

        outputs = self.expected_outputs(dataset)

        job = get_batch().new_job(self.name, self.get_job_attrs(dataset))
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.cpu(2).memory('highmem').storage('10Gi')
        job.command(f'annotate_dataset_small --input {str(input_mt)} --output {str(outputs)} ')
        return self.make_outputs(dataset, data=outputs, jobs=job)


@stage(required_stages=[AnnotateDatasetSmallVariantsWithHailQuery], analysis_type='custom', analysis_keys=['vcf'])
class AnnotatedDatasetMtToVcfWithHailQuery(DatasetStage):
    """
    Take the per-dataset annotated MT and write out as a VCF
    Optional stage set by dataset name in the config file
    """

    def expected_outputs(self, dataset: Dataset):
        """
        Expected to generate a VCF from the single-dataset MT
        """
        if family_sgs := get_family_sequencing_groups(dataset):
            return {
                'vcf': (
                    dataset.prefix()
                    / 'vcf'
                    / f'{get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.vcf.bgz'
                ),
                'index': (
                    dataset.prefix()
                    / 'vcf'
                    / f'{get_workflow().output_version}-{dataset.name}-{family_sgs["name_suffix"]}.vcf.bgz.tbi'
                ),
            }
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
        eligible_datasets = config_retrieve(['workflow', 'write_vcf'])
        if dataset.name not in eligible_datasets:
            return None

        mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetSmallVariantsWithHailQuery, key='mt')

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
        project_id=config_retrieve(['elasticsearch', 'password_project_id']),
        secret_name=config_retrieve(['elasticsearch', 'password_secret_id']),
        fail_gracefully=False,
    )


@stage(
    required_stages=[AnnotateDatasetSmallVariantsWithHailQuery],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'VARIANTS'},
)
class ExportMtAsEsIndex(DatasetStage):
    """
    Create a Seqr index.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        if family_sgs := get_family_sequencing_groups(dataset):
            index_name = (
                f'{dataset.name}-{sequencing_type}-{family_sgs["name_suffix"]}-{get_workflow().run_timestamp}'.lower()
            )
        else:
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
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetSmallVariantsWithHailQuery))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        outputs = self.expected_outputs(dataset)

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_bash_job(f'Generate {index_name} from {mt_path}')
        job._dirname = f'{index_name}-{job._token}'
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
