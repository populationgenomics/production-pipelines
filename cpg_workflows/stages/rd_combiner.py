from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.gcloud_composer import gcloud_compose_vcf_from_manifest
from cpg_workflows.jobs.rd_combiner.vqsr import (
    apply_recalibration_indels,
    apply_snp_vqsr_to_fragments,
    gather_tranches,
    train_vqsr_indels,
    train_vqsr_snp_tranches,
    train_vqsr_snps,
)
from cpg_workflows.targets import MultiCohort
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    MultiCohortStage,
    StageInput,
    StageOutput,
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


@stage(analysis_type='combiner', analysis_keys=['vds'])
class CreateVdsFromGvcfsWithHailCombiner(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path | str]:
        return {
            'vds': self.prefix / f'{multicohort.name}.vds',
            'combiner_plan': str(self.tmp_prefix / 'combiner_plan.json'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import combiner

        outputs: dict[str, str | Path] = self.expected_outputs(multicohort)

        # create these as empty lists instead of None, they have the same truthiness
        vds_path: str | None = None
        sg_ids_in_vds: list[str] = []

        if config_retrieve(['workflow', 'check_for_existing_vds'], True):
            # check for existing VDS
            get_logger(__file__).info('Checking for existing VDS')
            if existing_vds_analysis_entry := query_for_latest_vds(multicohort.analysis_dataset.name, 'combiner'):
                vds_path = existing_vds_analysis_entry['output']
                sg_ids_in_vds = [sg['id'] for sg in existing_vds_analysis_entry['sequencingGroups']]

        new_sg_gvcfs: list[str] = [
            str(sg.gvcf)
            for sg in multicohort.get_sequencing_groups()
            if (sg.gvcf is not None) and (sg.id not in sg_ids_in_vds)
        ]

        if not new_sg_gvcfs:
            get_logger(__file__).info('No GVCFs to combine')
            get_logger(__file__).info(f'Checking if VDS exists: {outputs["vds"]}: {outputs["vds"].exists()}')  # type: ignore
            return self.make_outputs(multicohort, outputs)

        j = get_batch().new_python_job('CreateVdsFromGvcfsWithHailCombiner', {'stage': self.name})
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.memory(config_retrieve(['combiner', 'memory']))
        j.storage(config_retrieve(['combiner', 'storage']))

        # Default to GRCh38 for reference if not specified
        j.call(
            combiner.run,
            output_vds_path=str(outputs['vds']),
            save_path=outputs['combiner_plan'],
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            tmp_prefix=str(self.tmp_prefix / 'temp_dir'),
            genome_build=genome_build(),
            gvcf_paths=new_sg_gvcfs,
            vds_paths=[vds_path] if vds_path else None,
        )

        return self.make_outputs(multicohort, outputs, j)


@stage(required_stages=[CreateVdsFromGvcfsWithHailCombiner], analysis_type='matrixtable', analysis_keys=['mt'])
class CreateDenseMtFromVdsWithHail(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        the MT and both shard_manifest files are Paths, so this stage will rerun if any of those are missing
        the VCFs are written as a directory, rather than a single VCF, so we can't check its existence well

        Needs a range of INFO fields to be present in the VCF
        """
        prefix = self.prefix

        return {
            'mt': prefix / f'{multicohort.name}.mt',
            # this will be the write path for fragments of sites-only VCF (header-per-shard)
            'hps_vcf_dir': str(prefix / f'{multicohort.name}.vcf.bgz'),
            # this will be the file which contains the name of all fragments (header-per-shard)
            'hps_shard_manifest': prefix / f'{multicohort.name}.vcf.bgz' / SHARD_MANIFEST,
            # this will be the write path for fragments of sites-only VCF (separate header)
            'separate_header_vcf_dir': str(prefix / f'{multicohort.name}_separate.vcf.bgz'),
            # this will be the file which contains the name of all fragments (separate header)
            'separate_header_manifest': prefix / f'{multicohort.name}_separate.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        output = self.expected_outputs(multicohort)

        j = get_batch().new_job('CreateDenseMtFromVdsWithHail')
        j.image(config_retrieve(['workflow', 'driver_image']))
        j.command(
            'mt_from_vds '
            f'--input {str(inputs.as_dict(multicohort, CreateVdsFromGvcfsWithHailCombiner)["vds"])} '
            f'--output {str(output["mt"])} '
            f'--sites_only {output["hps_vcf_dir"]} '
            f'--separate_header {output["separate_header_vcf_dir"]} ',
        )
        return self.make_outputs(multicohort, output, [j])


@stage(analysis_keys=['vcf'], analysis_type='vcf')
class ConcatenateVcfFragmentsWithGcloud(MultiCohortStage):
    """
    Takes a manifest of VCF fragments, and produces a single VCF file
    This is disconnected from the previous stage, but requires it to be run first
    So we check for the exact same output, and fail if we're not ready to start
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'vcf': self.prefix / 'gcloud_composed_sitesonly.vcf.bgz'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs to take a manifest of VCF fragments, and produce a single VCF file
        The VCF being composed here has a single header in a separate file, the first entry in the manifest
        This means we can combine the VCF header and data fragments through concatenation
        and the result will be a spec-compliant VCF
        """
        manifest_file = (
            multicohort.analysis_dataset.prefix()
            / 'rd_combiner'
            / get_workflow().output_version
            / 'CreateDenseMtFromVdsWithHail'
            / f'{multicohort.name}_separate.vcf.bgz'
            / SHARD_MANIFEST
        )

        if not manifest_file.exists():
            raise ValueError(f'Manifest file {str(manifest_file)} does not exist, run the rd_combiner workflow')

        outputs = self.expected_outputs(multicohort)

        jobs = gcloud_compose_vcf_from_manifest(
            manifest_path=manifest_file,
            intermediates_path=str(self.prefix / 'temporary_compose_intermediates'),
            output_path=str(outputs['vcf']),
            job_attrs={'stage': self.name},
        )

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(required_stages=[ConcatenateVcfFragmentsWithGcloud])
class TrainVqsrIndelModelOnCombinerData(MultiCohortStage):
    """
    Train VQSR Indel model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'indel_recalibrations': self.prefix / 'indel_recalibrations',
            'indel_tranches': self.prefix / 'indel_tranches',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_path(multicohort, ConcatenateVcfFragmentsWithGcloud, 'vcf')
        outputs = self.expected_outputs(multicohort)
        indel_calibration_job = train_vqsr_indels(
            sites_only_vcf=str(composed_sitesonly_vcf),
            indel_recal=str(outputs['indel_recalibrations']),
            indel_tranches=str(outputs['indel_tranches']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=indel_calibration_job)


@stage(required_stages=[ConcatenateVcfFragmentsWithGcloud])
class TrainVqsrSnpModelOnCombinerData(MultiCohortStage):
    """
    Train VQSR SNP model on the combiner data
    This is disconnected from the CreateDenseMtFromVdsWithHail stage, but requires it to be run first
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'snp_model': self.prefix / 'snp_model',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Submit jobs to train VQSR on the combiner data
        """

        composed_sitesonly_vcf = inputs.as_path(multicohort, ConcatenateVcfFragmentsWithGcloud, 'vcf')
        outputs = self.expected_outputs(multicohort)
        snp_calibration_job = train_vqsr_snps(
            sites_only_vcf=str(composed_sitesonly_vcf),
            snp_model=str(outputs['snp_model']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=snp_calibration_job)


@stage(required_stages=TrainVqsrSnpModelOnCombinerData)
class TrainVqsrSnpTranches(MultiCohortStage):
    """
    Scattered training of VQSR tranches for SNPs
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, str | Path]:

        return {
            'tranche_marker': self.tmp_prefix / 'tranches_trained',
            'temp_path': str(self.tmp_prefix / 'vqsr_snp_tranches'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        manifest_file = (
            multicohort.analysis_dataset.prefix()
            / 'rd_combiner'
            / get_workflow().output_version
            / 'CreateDenseMtFromVdsWithHail'
            / f'{multicohort.name}.vcf.bgz'
            / SHARD_MANIFEST
        )

        if not manifest_file.exists():
            raise ValueError(f'Manifest file {str(manifest_file)} does not exist, run the rd_combiner workflow')

        outputs = self.expected_outputs(multicohort)
        snp_model_path = inputs.as_path(target=multicohort, stage=TrainVqsrSnpModelOnCombinerData, key='snp_model')

        jobs = train_vqsr_snp_tranches(
            manifest_file=manifest_file,
            snp_model_path=str(snp_model_path),
            output_path=str(outputs['tranche_marker']),
            temp_path=to_path(outputs['temp_path']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(
    required_stages=[
        TrainVqsrSnpModelOnCombinerData,
        TrainVqsrSnpTranches,
    ],
)
class RunTrainedSnpVqsrOnCombinerFragments(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'vcf': self.prefix / 'vqsr.vcf.gz'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        manifest_file = (
            multicohort.analysis_dataset.prefix()
            / 'rd_combiner'
            / get_workflow().output_version
            / 'CreateDenseMtFromVdsWithHail'
            / f'{multicohort.name}.vcf.bgz'
            / SHARD_MANIFEST
        )

        if not manifest_file.exists():
            raise ValueError(f'Manifest file {manifest_file} does not exist, run the rd_combiner workflow')

        outputs = self.expected_outputs(multicohort)
        tranche_recal_temp = to_path(inputs.as_str(target=multicohort, stage=TrainVqsrSnpTranches, key='temp_path'))

        jobs = apply_snp_vqsr_to_fragments(
            manifest_file=manifest_file,
            temp_path=tranche_recal_temp,
            output_path=str(outputs['vcf']),
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(
    analysis_keys=['vcf'],
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

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'vcf': self.prefix / 'vqsr_snps_and_indels.vcf.gz'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        annotated_vcf = inputs.as_path(target=multicohort, stage=RunTrainedSnpVqsrOnCombinerFragments, key='vcf')
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
            output_path=outputs['vcf'],
            job_attrs={'stage': self.name},
        )
        return self.make_outputs(multicohort, data=outputs, jobs=indel_recal_job)
