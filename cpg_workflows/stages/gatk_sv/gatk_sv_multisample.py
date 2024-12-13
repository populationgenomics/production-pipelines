"""
All post-batching stages of the GATK-SV workflow
"""

from functools import cache
from itertools import combinations
from typing import Any

from google.api_core.exceptions import PermissionDenied

from cpg_utils import Path, to_path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, image_path, try_get_ar_guid
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch
from cpg_workflows.jobs import ploidy_table_from_ped
from cpg_workflows.jobs.gatk_sv import rename_sv_ids
from cpg_workflows.jobs.seqr_loader_sv import annotate_cohort_jobs_sv, annotate_dataset_jobs_sv
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    SV_CALLERS,
    CromwellJobSizes,
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_ref_panel,
    get_references,
    make_combined_ped,
    queue_annotate_strvctvre_job,
    queue_annotate_sv_jobs,
)
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.targets import Cohort, Dataset, MultiCohort
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    CohortStage,
    DatasetStage,
    MultiCohortStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

VCF_QUERY = gql(
    """
    query MyQuery($dataset: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: "sv"}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
            }
        }
    }
""",
)


@cache
def query_for_spicy_vcf(dataset: str) -> str | None:
    """
    query for the most recent previous SpiceUpSVIDs VCF
    the SpiceUpSVIDs Stage involves overwriting the generic sequential variant IDs
    with meaningful Identifiers, so we can track the same variant across different callsets

    Args:
        dataset (str): project to query for

    Returns:
        str, the path to the latest Spicy VCF
        or None, if there are no Spicy VCFs
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    result = query(VCF_QUERY, variables={'dataset': query_dataset})
    spice_by_date: dict[str, str] = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and analysis['output'].endswith('fresh_ids.vcf.bgz'):
            spice_by_date[analysis['timestampCompleted']] = analysis['output']

    if not spice_by_date:
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return spice_by_date[sorted(spice_by_date)[-1]]


@stage
class MakeCohortCombinedPed(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {'cohort_ped': self.get_stage_cohort_prefix(cohort) / 'ped_with_ref_panel.ped'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        output = self.expected_outputs(cohort)
        # write the pedigree, if it doesn't already exist
        if not to_path(output['cohort_ped']).exists():
            make_combined_ped(cohort, self.get_stage_cohort_prefix(cohort))

        return self.make_outputs(target=cohort, data=output)


def check_for_cohort_overlaps(multicohort: MultiCohort):
    """
    Check for overlapping cohorts in a MultiCohort.
    GATK-SV does not tolerate overlapping cohorts, so we check for this here.
    This is called once per MultiCohort, and raises an Exception if any overlaps are found

    Args:
        multicohort (MultiCohort): the MultiCohort to check
    """
    # placeholder for errors
    errors: list[str] = []
    sgs_per_cohort: dict[str, set[str]] = {}
    # grab all SG IDs per cohort
    for cohort in multicohort.get_cohorts():
        # shouldn't be possible, but guard against to be sure
        if cohort.id in sgs_per_cohort:
            raise ValueError(f'Cohort {cohort.id} already exists in {sgs_per_cohort}')

        # collect the SG IDs for this cohort
        sgs_per_cohort[cohort.id] = set(cohort.get_sequencing_group_ids())

    # pairwise iteration over cohort IDs
    for id1, id2 in combinations(sgs_per_cohort, 2):
        # if the IDs are the same, skip. Again, shouldn't be possible, but guard against to be sure
        if id1 == id2:
            continue
        # if there are overlapping SGs, raise an error
        if overlap := sgs_per_cohort[id1] & sgs_per_cohort[id2]:
            errors.append(f'Overlapping cohorts {id1} and {id2} have overlapping SGs: {overlap}')
    # upon findings any errors, raise an Exception and die
    if errors:
        raise ValueError('\n'.join(errors))


@stage
class MakeMultiCohortCombinedPed(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'multicohort_ped': self.prefix / 'ped_with_ref_panel.ped'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        # check that there are no overlapping cohorts
        check_for_cohort_overlaps(multicohort)

        output = self.expected_outputs(multicohort)
        # write the pedigree, if it doesn't already exist
        if not to_path(output['multicohort_ped']).exists():
            make_combined_ped(multicohort, self.prefix)

        return self.make_outputs(target=multicohort, data=output)


@stage(required_stages=[MakeCohortCombinedPed])
class GatherBatchEvidence(CohortStage):
    """
    This is the first Stage in the multisample GATK-SV workflow, running on a
    controlled subset of SGs as determined by the output of CreateSampleBatches.

    Using Analysis-Runner, include three additional config files:

    - configs/gatk_sv/stop_at_filter_batch.toml
        - this contains the instruction to stop prior to FilterBatch
    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - A custom config with the specific cohort in workflow.only_sgs

    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl

    it's critical to separate the ending with a dot, e.g.: `*.sr.txt.gz`,
    These files are passed to `gatk PrintSVEvidence`, that determines file
    format based on the file name.
    It would strongly expect the files to end exactly with either
    `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
    "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
    no suitable codecs found".
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """create the output paths for GatherBatchEvidence"""
        ending_by_key = {
            'cnmops_dup': 'DUP.header.bed.gz',
            'cnmops_dup_index': 'DUP.header.bed.gz.tbi',
            'cnmops_del': 'DEL.header.bed.gz',
            'cnmops_del_index': 'DEL.header.bed.gz.tbi',
            'cnmops_large_del': 'DEL.large.bed.gz',
            'cnmops_large_del_index': 'DEL.large.bed.gz.tb',
            'cnmops_large_dup': 'DUP.large.bed.gz',
            'cnmops_large_dup_index': 'DUP.large.bed.gz.tbi',
            'merged_SR': f'{self.name}.sr.txt.gz',
            'merged_SR_index': f'{self.name}.sr.txt.gz.tbi',
            'merged_PE': f'{self.name}.pe.txt.gz',
            'merged_PE_index': f'{self.name}.pe.txt.gz.tbi',
            'merged_BAF': f'{self.name}.baf.txt.gz',
            'merged_BAF_index': f'{self.name}.baf.txt.gz.tbi',
            'merged_bincov': f'{self.name}.RD.txt.gz',
            'merged_bincov_index': f'{self.name}.RD.txt.gz.tbi',
            'median_cov': 'medianCov.transposed.bed',
            'merged_dels': 'DEL.bed.gz',
            'merged_dups': 'DUP.bed.gz',
        }

        # we don't run metrics as standard, only expect the output if we choose to run
        if config_retrieve(['resource_overrides', self.name, 'run_matrix_qc'], False):
            ending_by_key.update(
                {
                    'Matrix_QC_plot': '00_matrix_FC_QC.png',
                    'SR_stats': 'SR.QC_matrix.txt',
                    'PE_stats': 'PE.QC_matrix.txt',
                    'BAF_stats': 'BAF.QC_matrix.txt',
                    'RD_stats': 'RD.QC_matrix.txt',
                },
            )

        for caller in SV_CALLERS:
            ending_by_key[f'std_{caller}_vcf_tar'] = f'{caller}.tar.gz'

        return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """Add jobs to Batch"""
        sequencing_groups = cohort.get_sequencing_groups(only_active=True)
        pedigree_input = inputs.as_path(target=cohort, stage=MakeCohortCombinedPed, key='cohort_ped')

        input_dict: dict[str, Any] = {
            'batch': cohort.id,
            'samples': [sg.id for sg in sequencing_groups],
            'ped_file': str(pedigree_input),
            'counts': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.coverage_counts.tsv.gz')
                for sequencing_group in sequencing_groups
            ],
            'SR_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sr.txt.gz')
                for sequencing_group in sequencing_groups
            ],
            'PE_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.pe.txt.gz')
                for sequencing_group in sequencing_groups
            ],
            'SD_files': [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sd.txt.gz')
                for sequencing_group in sequencing_groups
            ],
            'ref_copy_number_autosomal_contigs': 2,
            'allosomal_contigs': ['chrX', 'chrY'],
            'gcnv_qs_cutoff': 30,
            'min_svsize': 50,
            'run_matrix_qc': True,
            'matrix_qc_distance': 1000000,
            'ref_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.{caller}.vcf.gz')
                for sequencing_group in sequencing_groups
            ]

        input_dict |= get_references(
            [
                'genome_file',
                'primary_contigs_fai',
                {'sd_locs_vcf': 'dbsnp_vcf'},
                {'cnmops_chrom_file': 'autosome_file'},
                'cnmops_exclude_list',
                {'cnmops_allo_file': 'allosome_file'},
                'cytoband',
                'mei_bed',
            ],
        )

        # reference panel gCNV models
        input_dict |= get_ref_panel()

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
                'cnmops_docker',
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # this step runs for approximately 15 hours
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.LARGE,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeCohortCombinedPed, GatherBatchEvidence])
class ClusterBatch(CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#clusterbatch
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        * Clustered SV VCFs
        * Clustered depth-only call VCF
        """

        ending_by_key = {}

        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'clustered_{caller}_vcf'] = f'clustered-{caller}.vcf.gz'
            ending_by_key[f'clustered_{caller}_vcf_index'] = f'clustered-{caller}.vcf.gz.tbi'
        return {key: self.get_stage_cohort_prefix(cohort) / fname for key, fname in ending_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """
        batch_evidence_d = inputs.as_dict(cohort, GatherBatchEvidence)
        pedigree_input = inputs.as_path(target=cohort, stage=MakeCohortCombinedPed, key='cohort_ped')

        input_dict: dict[str, Any] = {
            'batch': cohort.id,
            'del_bed': str(batch_evidence_d['merged_dels']),
            'dup_bed': str(batch_evidence_d['merged_dups']),
            'ped_file': str(pedigree_input),
            'depth_exclude_overlap_fraction': 0.5,
            'depth_interval_overlap': 0.8,
            'depth_clustering_algorithm': 'SINGLE_LINKAGE',
            'pesr_interval_overlap': 0.1,
            'pesr_breakend_window': 300,
            'pesr_clustering_algorithm': 'SINGLE_LINKAGE',
            'reference_fasta': str(get_fasta()),
            'reference_fasta_fai': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcf_tar'] = str(batch_evidence_d[f'std_{caller}_vcf_tar'])

        input_dict |= get_images(
            [
                'linux_docker',
                'sv_pipeline_base_docker',
                'gatk_docker',
                'sv_base_mini_docker',
                'sv_pipeline_docker',
            ],
        )

        input_dict |= get_references(
            [
                {'contig_list': 'primary_contigs_list'},
                {'depth_exclude_intervals': 'depth_exclude_list'},
                {'pesr_exclude_intervals': 'pesr_exclude_list'},
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for approx 1 hour
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeCohortCombinedPed, ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(CohortStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Metrics files
        """

        return {
            'metrics': self.get_stage_cohort_prefix(cohort) / 'metrics.tsv',
            'metrics_common': self.get_stage_cohort_prefix(cohort) / 'metrics_common.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)
        gatherbatchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)
        pedigree_input = inputs.as_path(target=cohort, stage=MakeCohortCombinedPed, key='cohort_ped')

        input_dict: dict[str, Any] = {
            'batch': cohort.id,
            'baf_metrics': gatherbatchevidence_d['merged_BAF'],
            'discfile': gatherbatchevidence_d['merged_PE'],
            'coveragefile': gatherbatchevidence_d['merged_bincov'],
            'splitfile': gatherbatchevidence_d['merged_SR'],
            'medianfile': gatherbatchevidence_d['median_cov'],
            'BAF_split_size': 10000,
            'RD_split_size': 10000,
            'PE_split_size': 10000,
            'SR_split_size': 1000,
            'common_cnv_size_cutoff': 5000,
            'ped_file': str(pedigree_input),
            'ref_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_pipeline_rdtest_docker',
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_base_docker',
                'linux_docker',
            ],
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
                'rmsk',
                'segdups',
                {'autosome_contigs': 'autosome_file'},
                {'allosome_contigs': 'allosome_file'},
            ],
        )

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for approx 4-5 hours
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeCohortCombinedPed, GenerateBatchMetrics, ClusterBatch])
class FilterBatch(CohortStage):
    """
    Filters poor quality variants and filters outlier samples.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        * Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        * Filtered depth-only call VCF with outlier samples excluded
        * Random forest cutoffs file
        * PED file with outlier samples excluded
        """

        ending_by_key: dict = {
            'filtered_pesr_vcf': 'filtered_pesr_merged.vcf.gz',
            'cutoffs': 'cutoffs',
            'scores': 'updated_scores',
            'RF_intermediate_files': 'RF_intermediate_files.tar.gz',
            'outlier_samples_excluded_file': 'outliers.samples.list',
            'batch_samples_postOutlierExclusion_file': 'outliers_excluded.samples.list',
        }

        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'filtered_{caller}_vcf'] = f'filtered-{caller}.vcf.gz'

            # unsure why, scramble doesn't export this file
            if caller != 'scramble':
                ending_by_key[f'sites_filtered_{caller}_vcf'] = f'sites-filtered-{caller}.vcf.gz'

        ending_by_key['sv_counts'] = [f'{caller}.with_evidence.svcounts.txt' for caller in SV_CALLERS + ['depth']]
        ending_by_key['sv_count_plots'] = [
            f'{caller}.with_evidence.all_SVTYPEs.counts_per_sample.png' for caller in SV_CALLERS + ['depth']
        ]
        d: dict[str, Path | list[Path]] = {}
        for key, ending in ending_by_key.items():
            if isinstance(ending, str):
                d[key] = self.get_stage_cohort_prefix(cohort) / ending
            elif isinstance(ending, list):
                d[key] = [self.get_stage_cohort_prefix(cohort) / e for e in ending]
        return d

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        metrics_d = inputs.as_dict(cohort, GenerateBatchMetrics)
        clusterbatch_d = inputs.as_dict(cohort, ClusterBatch)
        pedigree_input = inputs.as_path(target=cohort, stage=MakeCohortCombinedPed, key='cohort_ped')

        input_dict: dict[str, Any] = {
            'batch': cohort.id,
            'ped_file': str(pedigree_input),
            'evidence_metrics': metrics_d['metrics'],
            'evidence_metrics_common': metrics_d['metrics_common'],
            'outlier_cutoff_nIQR': '6',
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            ['sv_pipeline_base_docker', 'sv_pipeline_docker', 'sv_base_mini_docker', 'linux_docker'],
        )

        input_dict |= get_references(['primary_contigs_list'])

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for over an hour
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=FilterBatch)
class MergeBatchSites(MultiCohortStage):
    """
    OMG, my first true MultiCohortStage!

    This Stage runs between GenerateBatchMetrics and FilterBatch
    This takes the component VCFs from individual batches and merges them into
    a single VCF (one for PESR and one for depth-only calls).

    This only has to be run once for all Cohorts, and technically doesn't
    need to be associated with a group of SGs. However, including one (via the
    config containing `only_sgs`) allows us to keep the output VCFs within the
    output structure generated by previous stages in this run.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        generate them there outputs
        """
        return {
            'cohort_pesr_vcf': self.prefix / 'cohort_pesr.vcf.gz',
            'cohort_depth_vcf': self.prefix / 'cohort_depth.vcf.gz',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        generate a MergeBatchSites job
        """

        # take from previous per-cohort outputs
        filter_batch_outputs = inputs.as_dict_by_target(FilterBatch)
        pesr_vcfs = [filter_batch_outputs[cohort.id]['filtered_pesr_vcf'] for cohort in multicohort.get_cohorts()]
        depth_vcfs = [filter_batch_outputs[cohort.id]['filtered_depth_vcf'] for cohort in multicohort.get_cohorts()]

        input_dict: dict = {'cohort': multicohort.name, 'depth_vcfs': depth_vcfs, 'pesr_vcfs': pesr_vcfs}
        input_dict |= get_images(['sv_pipeline_docker'])
        expected_d = self.expected_outputs(multicohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(analysis_type='sv', analysis_keys=['exclusion_list'], required_stages=FilterBatch)
class CombineExclusionLists(MultiCohortStage):
    """
    Takes the per-batch lists of excluded sample IDs and combines
    them into a single file for use in the SV pipeline

    This will be used to remove any filtered samples from consideration in
    subsequent stages, and to remove the CPG ID registration
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Create dictionary of names -> output paths
        This variable is a Path to make sure it gets existence checked
        We need this quick stage to run each time
        """

        return {'exclusion_list': self.prefix / 'combined_exclusion_list.txt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        queue job to combine exclusion lists
        """

        filter_batch_outputs = inputs.as_dict_by_target(FilterBatch)
        all_filter_lists = [
            str(filter_batch_outputs[cohort.id]['outlier_samples_excluded_file'])
            for cohort in multicohort.get_cohorts()
        ]

        output = self.expected_outputs(multicohort)
        job = get_batch().new_job('Concatenate all sample exclusion files')
        job.image(config_retrieve(key=['workflow', 'driver_image']))
        authenticate_cloud_credentials_in_job(job)
        job.command(f'gcloud storage objects compose {" ".join(all_filter_lists)} {output["exclusion_list"]}')
        return self.make_outputs(multicohort, data=output, jobs=job)


@stage(
    required_stages=[FilterBatch, GatherBatchEvidence, MergeBatchSites],
    analysis_type='sv',
    analysis_keys=['genotyped_depth_vcf', 'genotyped_pesr_vcf'],
)
class GenotypeBatch(CohortStage):
    """
    The final Cohort Stage - executed for each individual Batch
    Genotypes a batch of samples across filtered variants combined across all batches.
    Include the additional config file:
    - configs/gatk_sv/use_for_all_workflows.toml; contains all required images and references
    The final CohortStage in this workflow
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """

        ending_by_key = {
            'sr_bothside_pass': 'genotype_SR_part2_bothside_pass.txt',
            'sr_background_fail': 'genotype_SR_part2_background_fail.txt',
            'trained_PE_metrics': 'pe_metric_file.txt',
            'trained_SR_metrics': 'sr_metric_file.txt',
            'regeno_coverage_medians': 'regeno.coverage_medians_merged.bed',
        }

        for mode in ['pesr', 'depth']:
            ending_by_key |= {
                f'trained_genotype_{mode}_pesr_sepcutoff': f'{mode}.pesr_sepcutoff.txt',
                f'trained_genotype_{mode}_depth_sepcutoff': f'{mode}.depth_sepcutoff.txt',
                f'genotyped_{mode}_vcf': f'{mode}.vcf.gz',
                f'genotyped_{mode}_vcf_index': f'{mode}.vcf.gz.tbi',
            }

        return {
            key: self.get_stage_cohort_prefix(cohort) / get_workflow().output_version / fname
            for key, fname in ending_by_key.items()
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        filterbatch_d = inputs.as_dict(cohort, FilterBatch)
        batchevidence_d = inputs.as_dict(cohort, GatherBatchEvidence)

        # workaround for mypy - cohort.multicohort is MultiCohort | None, and we require not-None for the inputs call
        this_multicohort = cohort.multicohort
        assert this_multicohort is not None, 'Multicohort cannot be None'
        mergebatch_d = inputs.as_dict(this_multicohort, MergeBatchSites)

        input_dict: dict[str, Any] = {
            'batch': cohort.id,
            'n_per_split': config_retrieve(['resource_overrides', 'GenotypeBatch', 'n_per_split'], 5000),
            'n_RD_genotype_bins': config_retrieve(
                ['resource_overrides', 'GenotypeBatch', 'n_RD_genotype_bins'],
                100000,
            ),
            'coveragefile': batchevidence_d['merged_bincov'],
            'coveragefile_index': batchevidence_d['merged_bincov_index'],
            'discfile': batchevidence_d['merged_PE'],
            'discfile_index': batchevidence_d['merged_PE_index'],
            'splitfile': batchevidence_d['merged_SR'],
            'splitfile_index': batchevidence_d['merged_SR_index'],
            'medianfile': batchevidence_d['median_cov'],
            'rf_cutoffs': filterbatch_d['cutoffs'],
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'reference_build': 'hg38',
            'batch_depth_vcf': filterbatch_d['filtered_depth_vcf'],
            'batch_pesr_vcf': filterbatch_d['filtered_pesr_vcf'],
            'cohort_depth_vcf': mergebatch_d['cohort_depth_vcf'],
            'cohort_pesr_vcf': mergebatch_d['cohort_pesr_vcf'],
        }

        input_dict |= get_images(
            [
                'sv_pipeline_base_docker',
                'sv_base_mini_docker',
                'sv_pipeline_docker',
                'sv_pipeline_rdtest_docker',
                'linux_docker',
            ],
        )
        input_dict |= get_references(['primary_contigs_list', 'bin_exclude', 'seed_cutoffs', 'pesr_exclude_list'])
        # 2 additional (optional) references CPG doesn't have:
        # - sr_hom_cutoff_multiplier
        # - sr_median_hom_ins

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeMultiCohortCombinedPed, GatherBatchEvidence, GenotypeBatch, FilterBatch])
class MakeCohortVcf(MultiCohortStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """create output paths"""
        return {
            'vcf': self.prefix / 'cleaned.vcf.gz',
            'vcf_index': self.prefix / 'cleaned.vcf.gz.tbi',
            'vcf_qc': self.prefix / 'cleaned_SV_VCF_QC_output.tar.gz',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Instead of taking a direct dependency on the previous stage,
        we use the output hash to find all the previous batches
        """
        gatherbatchevidence_outputs = inputs.as_dict_by_target(GatherBatchEvidence)
        genotypebatch_outputs = inputs.as_dict_by_target(GenotypeBatch)
        filterbatch_outputs = inputs.as_dict_by_target(FilterBatch)
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')

        # get the names of all contained cohorts
        all_batch_names: list[str] = [cohort.id for cohort in multicohort.get_cohorts()]

        pesr_vcfs = [genotypebatch_outputs[cohort]['genotyped_pesr_vcf'] for cohort in all_batch_names]
        depth_vcfs = [genotypebatch_outputs[cohort]['genotyped_depth_vcf'] for cohort in all_batch_names]
        sr_pass = [genotypebatch_outputs[cohort]['sr_bothside_pass'] for cohort in all_batch_names]
        sr_fail = [genotypebatch_outputs[cohort]['sr_background_fail'] for cohort in all_batch_names]
        depth_depth_cutoff = [
            genotypebatch_outputs[cohort]['trained_genotype_depth_depth_sepcutoff'] for cohort in all_batch_names
        ]
        filter_batch_cutoffs = [filterbatch_outputs[cohort]['cutoffs'] for cohort in all_batch_names]
        bincov_files = [gatherbatchevidence_outputs[cohort]['merged_bincov'] for cohort in all_batch_names]
        disc_files = [gatherbatchevidence_outputs[cohort]['merged_PE'] for cohort in all_batch_names]
        median_cov_files = [gatherbatchevidence_outputs[cohort]['median_cov'] for cohort in all_batch_names]

        input_dict: dict[str, Any] = {
            'cohort_name': multicohort.name,
            'batches': all_batch_names,
            'ped_file': str(pedigree_input),
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'chr_x': 'chrX',
            'chr_y': 'chrY',
            'min_sr_background_fail_batches': 0.5,
            'max_shard_size_resolve': 500,
            'max_shards_per_chrom_clean_vcf_step1': 200,
            'min_records_per_shard_clean_vcf_step1': 5000,
            'clean_vcf1b_records_per_shard': 10000,
            'samples_per_clean_vcf_step2_shard': 100,
            'clean_vcf5_records_per_shard': 5000,
            'random_seed': 0,
            # not explicit, but these VCFs require indices
            'pesr_vcfs': pesr_vcfs,
            'depth_vcfs': depth_vcfs,
            'disc_files': disc_files,
            'bincov_files': bincov_files,
            'raw_sr_bothside_pass_files': sr_pass,
            'raw_sr_background_fail_files': sr_fail,
            'depth_gt_rd_sep_files': depth_depth_cutoff,
            'median_coverage_files': median_cov_files,
            'rf_cutoff_files': filter_batch_cutoffs,
        }

        input_dict |= get_references(
            [
                'bin_exclude',
                'mei_bed',
                'depth_exclude_list',
                'empty_file',
                # same attr, two names
                'primary_contigs_list',
                {'contig_list': 'primary_contigs_list'},
                {'allosome_fai': 'allosome_file'},
                {'cytobands': 'cytoband'},
                {'pe_exclude_list': 'pesr_exclude_list'},
            ],
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_pipeline_base_docker',
                'sv_pipeline_hail_docker',
                'sv_pipeline_updates_docker',
                'sv_pipeline_rdtest_docker',
                'sv_pipeline_qc_docker',
                'sv_base_mini_docker',
                'linux_docker',
            ],
        )
        expected_d = self.expected_outputs(multicohort)

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
            job_size=CromwellJobSizes.MEDIUM,
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeMultiCohortCombinedPed, MakeCohortVcf])
class FormatVcfForGatk(MultiCohortStage):

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
        input_dict: dict[str, Any] = {
            'prefix': multicohort.name,
            'vcf': inputs.as_dict(multicohort, MakeCohortVcf)['vcf'],
            'ped_file': str(pedigree_input),
        }
        input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(multicohort)

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
            job_size=CromwellJobSizes.SMALL,
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeMultiCohortCombinedPed, ClusterBatch, MakeCohortVcf])
class JoinRawCalls(MultiCohortStage):
    """
    Joins all individually clustered caller results
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'joined_raw_calls_vcf': self.prefix / 'raw_clustered_calls.vcf.gz',
            'joined_raw_calls_vcf_index': self.prefix / 'raw_clustered_calls.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
        input_dict: dict[str, Any] = {
            'FormatVcfForGatk.formatter_args': '--fix-end',
            'prefix': multicohort.name,
            'ped_file': str(pedigree_input),
            'reference_fasta': get_fasta(),
            'reference_fasta_fai': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
        }
        input_dict |= get_images(['gatk_docker', 'sv_pipeline_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        # add all clustered _caller_ files, plus indices
        clusterbatch_outputs = inputs.as_dict_by_target(ClusterBatch)

        # get the names of all contained cohorts
        all_batch_names: list[str] = [cohort.id for cohort in multicohort.get_cohorts()]
        for caller in SV_CALLERS + ['depth']:
            input_dict[f'clustered_{caller}_vcfs'] = [
                clusterbatch_outputs[cohort][f'clustered_{caller}_vcf'] for cohort in all_batch_names
            ]
            input_dict[f'clustered_{caller}_vcf_indexes'] = [
                clusterbatch_outputs[cohort][f'clustered_{caller}_vcf_index'] for cohort in all_batch_names
            ]

        expected_d = self.expected_outputs(multicohort)

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(required_stages=[JoinRawCalls, FormatVcfForGatk])
class SVConcordance(MultiCohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'concordance_vcf': self.prefix / 'sv_concordance.vcf.gz',
            'concordance_vcf_index': self.prefix / 'sv_concordance.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV concordance
        """

        input_dict: dict[str, Any] = {
            'output_prefix': multicohort.name,
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'eval_vcf': inputs.as_dict(multicohort, FormatVcfForGatk)['gatk_formatted_vcf'],
            'truth_vcf': inputs.as_dict(multicohort, JoinRawCalls)['joined_raw_calls_vcf'],
        }
        input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(multicohort)

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(required_stages=[MakeMultiCohortCombinedPed, SVConcordance])
class GeneratePloidyTable(MultiCohortStage):
    """
    Quick PythonJob to generate a ploidy table
    Calls a homebrewed version of this table generator:
    github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/ploidy_table_from_ped.py
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        only one output, the ploidy table
        """

        return {'ploidy_table': self.prefix / 'ploidy_table.txt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')

        py_job = get_batch().new_python_job('create_ploidy_table')
        py_job.image(config_retrieve(['workflow', 'driver_image']))

        contig_path = get_references(['primary_contigs_list'])['primary_contigs_list']

        expected_d = self.expected_outputs(multicohort)
        py_job.call(
            ploidy_table_from_ped.generate_ploidy_table,
            str(pedigree_input),
            contig_path,
            str(expected_d['ploidy_table']),
        )

        return self.make_outputs(multicohort, data=expected_d, jobs=py_job)


@stage(required_stages=[MakeMultiCohortCombinedPed, GeneratePloidyTable, SVConcordance])
class FilterGenotypes(MultiCohortStage):
    """
    Steps required to post-filter called genotypes
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'filtered_vcf': self.prefix / 'filtered.vcf.gz',
            'filtered_vcf_index': self.prefix / 'filtered.vcf.gz.tbi',
            'unfiltered_recalibrated_vcf': self.prefix / 'unfiltered_recalibrated.vcf.gz',
            'unfiltered_recalibrated_vcf_index': self.prefix / 'unfiltered_recalibrated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        pedigree_input = inputs.as_path(target=multicohort, stage=MakeMultiCohortCombinedPed, key='multicohort_ped')
        input_dict = {
            'output_prefix': multicohort.name,
            'vcf': inputs.as_dict(multicohort, SVConcordance)['concordance_vcf'],
            'ploidy_table': inputs.as_dict(multicohort, GeneratePloidyTable)['ploidy_table'],
            'ped_file': str(pedigree_input),
            'fmax_beta': config_retrieve(['references', 'gatk_sv', 'fmax_beta'], 0.3),
            'recalibrate_gq_args': config_retrieve(['references', 'gatk_sv', 'recalibrate_gq_args']),
            'sl_filter_args': config_retrieve(['references', 'gatk_sv', 'sl_filter_args']),
        }
        assert input_dict['recalibrate_gq_args'] and input_dict['sl_filter_args'], input_dict

        input_dict |= get_images(['linux_docker', 'sv_base_mini_docker', 'sv_pipeline_docker'])
        # use a non-standard GATK image containing required filtering tool
        input_dict['gatk_docker'] = get_images(['gq_recalibrator_docker'])['gq_recalibrator_docker']
        input_dict |= get_references(
            [{'gq_recalibrator_model_file': 'aou_filtering_model'}, 'primary_contigs_fai'],
        )

        # something a little trickier - we need to get various genome tracks
        input_dict['genome_tracks'] = list(
            get_references(config_retrieve(['references', 'gatk_sv', 'genome_tracks'], [])).values(),
        )

        expected_d = self.expected_outputs(multicohort)

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(required_stages=FilterGenotypes)
class UpdateStructuralVariantIDs(MultiCohortStage):
    """
    Runs SVConcordance between the results of this callset and the results of a previous callset
    This causes the Variant IDs of a matching variant to be updated to the previous callset's ID
    Consistency of Variant ID is crucial to Seqr/AIP identifying the same variant across different callsets
    By default GATK-SV creates an auto-incrementing ID, rather than one based on variant attributes
    If a new call is added at a chromosome start for any sample, the ID for all variants on the chromosome
    will be shifted up, meaning that Seqr labels/any other process linked to the ID will be incorrect
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'concordance_vcf': self.prefix / 'updated_ids.vcf.gz',
            'concordance_vcf_index': self.prefix / 'updated_ids.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:

        # allow for no prior name/IDs
        if not (spicy_vcf := query_for_spicy_vcf(multicohort.analysis_dataset.name)):
            get_logger().info('No previous Spicy VCF found for {cohort.analysis_dataset.name}')
            return self.make_outputs(multicohort, skipped=True)

        expected_d = self.expected_outputs(multicohort)

        # run concordance between this and the previous VCF
        input_dict: dict = {
            'output_prefix': multicohort.name,
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'eval_vcf': inputs.as_path(multicohort, FilterGenotypes, key='filtered_vcf'),
            'truth_vcf': spicy_vcf,
        }
        input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        jobs = add_gatk_sv_jobs(
            dataset=multicohort.analysis_dataset,
            wfl_name='SVConcordance',
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels={'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()},
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=jobs)


@stage(
    required_stages=[FilterGenotypes, UpdateStructuralVariantIDs],
    analysis_type='sv',
    analysis_keys=['wham_filtered_vcf'],
)
class FilterWham(MultiCohortStage):
    """
    Filters the VCF to remove deletions only called by Wham
    github.com/broadinstitute/gatk-sv/blob/main/wdl/ApplyManualVariantFilter.wdl
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        # index here is implicit
        return {'wham_filtered_vcf': self.prefix / 'filtered.vcf.bgz'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        # read the concordance-with-prev-batch VCF if appropriate, otherwise use the filtered VCF
        if query_for_spicy_vcf(multicohort.analysis_dataset.name):
            get_logger().info(f'Variant IDs were updated for {multicohort.analysis_dataset.name}')
            input_vcf = inputs.as_dict(multicohort, UpdateStructuralVariantIDs)['concordance_vcf']
        else:
            get_logger().info(f'No Spicy VCF was found, default IDs for {multicohort.analysis_dataset.name}')
            input_vcf = inputs.as_dict(multicohort, FilterGenotypes)['filtered_vcf']

        in_vcf = get_batch().read_input_group(**{'vcf.gz': input_vcf, 'vcf.gz.tbi': f'{input_vcf}.tbi'})['vcf.gz']
        job = get_batch().new_job('Filter Wham', attributes={'tool': 'bcftools'})
        job.image(image_path('bcftools')).cpu(1).memory('highmem').storage('20Gi')
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        job.command(
            'bcftools view -e \'SVTYPE=="DEL" && COUNT(ALGORITHMS)==1 && ALGORITHMS=="wham"\' '
            f'{in_vcf} | bgzip -c  > {job.output["vcf.bgz"]}',
        )
        job.command(f'tabix {job.output["vcf.bgz"]}')
        get_batch().write_output(
            job.output,
            str(self.expected_outputs(multicohort)['wham_filtered_vcf']).replace('.vcf.bgz', ''),
        )

        expected_d = self.expected_outputs(multicohort)
        return self.make_outputs(multicohort, data=expected_d, jobs=job)


@stage(required_stages=[FilterWham], analysis_type='sv', analysis_keys=['annotated_vcf'])
class AnnotateVcf(MultiCohortStage):
    """
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.

    Note: the annotation stage is stupid, and re-orders the VCF by ID instead of position
    This means that we run SVConcordance before this stage, to ensure that the IDs are correct
    But we don't apply those IDs, in case multiple variants map to the same ID, which would
    cause the variants to become unsorted when ordered by ID.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        return {
            'annotated_vcf': self.prefix / 'filtered_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'filtered_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Configure and queue jobs for SV annotation. Passing the VCF Index has become implicit
        """
        input_vcf = inputs.as_dict(multicohort, FilterWham)['wham_filtered_vcf']
        expected_out = self.expected_outputs(multicohort)
        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}
        job_or_none = queue_annotate_sv_jobs(multicohort, self.prefix, input_vcf, expected_out, billing_labels)
        return self.make_outputs(multicohort, data=expected_out, jobs=job_or_none)


@stage(required_stages=AnnotateVcf)
class AnnotateVcfWithStrvctvre(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.prefix / 'strvctvre_annotated.vcf.gz',
            'strvctvre_vcf_index': self.prefix / 'strvctvre_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        input_dict = inputs.as_dict(multicohort, AnnotateVcf)
        expected_d = self.expected_outputs(multicohort)

        # read vcf and index into the batch
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']

        strvctvre_job = queue_annotate_strvctvre_job(input_vcf, str(expected_d['strvctvre_vcf']), self.get_job_attrs())

        return self.make_outputs(multicohort, data=expected_d, jobs=strvctvre_job)


@stage(required_stages=AnnotateVcfWithStrvctvre, analysis_type='sv', analysis_keys=['new_id_vcf'])
class SpiceUpSVIDs(MultiCohortStage):
    """
    Overwrites the GATK-SV assigned IDs with a meaningful ID
    This new ID is either taken from an equivalent variant ID in the previous callset (found through SVConcordance)
    or a new one is generated based on the call attributes itself
    A boolean flag is used to switch behaviour, based on the truthiness of the return value of query_for_spicy_vcf
    If a VCF is returned, True, which also means a VCF should have been used in UpdateStructuralVariantIDs
    If the return is None, False, and the IDs will be generated from the call attributes alone
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'new_id_vcf': self.prefix / 'fresh_ids.vcf.bgz',
            'new_id_index': self.prefix / 'fresh_ids.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        input_vcf = get_batch().read_input(str(inputs.as_path(multicohort, AnnotateVcfWithStrvctvre, 'strvctvre_vcf')))
        expected_output = self.expected_outputs(multicohort)

        # update the IDs using a PythonJob
        pyjob = get_batch().new_python_job('rename_sv_ids')
        pyjob.storage('10Gi')
        skip_prior_names = bool(query_for_spicy_vcf(multicohort.analysis_dataset.name))
        pyjob.call(rename_sv_ids, input_vcf, pyjob.output, skip_prior_names)

        # then compress & run tabix on that plain text result
        bcftools_job = get_batch().new_job('bgzip and tabix')
        bcftools_job.image(image_path('bcftools'))
        bcftools_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        bcftools_job.command(f'bcftools view {pyjob.output} | bgzip -c > {bcftools_job.output["vcf.bgz"]}')
        bcftools_job.command(f'tabix {bcftools_job.output["vcf.bgz"]}')  # type: ignore

        # get the output root to write to
        get_batch().write_output(bcftools_job.output, str(expected_output['new_id_vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(multicohort, data=expected_output, jobs=[pyjob, bcftools_job])


@stage(required_stages=SpiceUpSVIDs)
class AnnotateCohortSv(MultiCohortStage):
    """
    First step to transform annotated SV callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        Expected to write a matrix table.
        """
        return {'tmp_prefix': str(self.tmp_prefix), 'mt': self.prefix / 'cohort_sv.mt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)

        vcf_path = inputs.as_path(target=multicohort, stage=SpiceUpSVIDs, key='new_id_vcf')
        checkpoint_prefix = to_path(outputs['tmp_prefix']) / 'checkpoints'

        job = annotate_cohort_jobs_sv(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            checkpoint_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[CombineExclusionLists, AnnotateCohortSv], analysis_type='sv', analysis_keys=['mt'])
class AnnotateDatasetSv(DatasetStage):
    """
    Subset the MT to be this Dataset only
    Then work up all the genotype values
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to generate a matrix table
        """
        return {'mt': (dataset.prefix() / 'mt' / f'SV-{get_workflow().output_version}-{dataset.name}.mt')}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        # only create dataset MTs for datasets specified in the config
        eligible_datasets = config_retrieve(['workflow', 'write_mt_for_datasets'], default=[])
        if dataset.name not in eligible_datasets:
            get_logger().info(f'Skipping AnnotateDatasetSv mt subsetting for {dataset}')
            return None

        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortSv, key='mt')
        exclusion_file = inputs.as_path(get_multicohort(), stage=CombineExclusionLists, key='exclusion_list')

        outputs = self.expected_outputs(dataset)

        jobs = annotate_dataset_jobs_sv(
            mt_path=mt_path,
            sgids=dataset.get_sequencing_group_ids(),
            out_mt_path=outputs['mt'],
            tmp_prefix=self.tmp_prefix / dataset.name / 'checkpoints',
            job_attrs=self.get_job_attrs(dataset),
            exclusion_file=str(exclusion_file),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage(
    required_stages=[AnnotateDatasetSv],
    analysis_type='es-index',  # specific type of es index
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SV'},
)
class MtToEsSv(DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-SV-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses the non-DataProc MT-to-ES conversion script
        """
        # only create the elasticsearch index for the datasets specified in the config
        eligible_datasets = config_retrieve(['workflow', 'create_es_index_for_datasets'], default=[])
        if dataset.name not in eligible_datasets:
            get_logger().info(f'Skipping SV ES index creation for {dataset}')
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

        outputs = self.expected_outputs(dataset)

        # get the absolute path to the MT
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetSv, key='mt'))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_job(f'Generate {index_name} from {mt_path}')
        # set all job attributes in one bash
        job.cpu(4).memory('lowmem').storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))

        # localise the MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

        # run the export from the localised MT - this job writes no new data, just transforms and exports over network
        job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
