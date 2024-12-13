"""
Stages that implement GATK-gCNV.
"""

import json
from functools import lru_cache

from google.api_core.exceptions import PermissionDenied

from cpg_utils import Path, to_path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, image_path, reference_path, try_get_ar_guid
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import gcnv
from cpg_workflows.stages.gatk_sv.gatk_sv_common import get_images, get_references, queue_annotate_sv_jobs
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.targets import Cohort, Dataset, MultiCohort, SequencingGroup
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    CohortStage,
    DatasetStage,
    MultiCohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)


@lru_cache(maxsize=None)
def get_cohort_for_sgid(sgid: str) -> Cohort:
    """
    Return the cohort that contains this sgid
    until a better central method exists this needs to run multiple times
    so build a method here with caching
    """
    for c in get_multicohort().get_cohorts():
        if sgid in c.get_sequencing_group_ids():
            return c
    raise ValueError(f'Could not find cohort for {sgid}')


@stage
class SetSGIDOrdering(CohortStage):
    """
    Set the order of the sequencing groups in the cohort
    Push this to a file _now_, and use it later
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {'sgid_order': self.get_stage_cohort_prefix(cohort) / 'sgid_order.json'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        sorted_sgids = sorted(cohort.get_sequencing_group_ids())
        with self.expected_outputs(cohort)['sgid_order'].open('w') as f_handler:
            json.dump(sorted_sgids, f_handler, indent=2)
        return self.make_outputs(cohort, data=self.expected_outputs(cohort))


@stage
class PrepareIntervals(MultiCohortStage):
    """
    Interval preparation steps that don't require the sample read counts:
    PreprocessIntervals and AnnotateIntervals.
    This is a multicohort stage - we only ever co-process SGIDs on a unified capture
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        prefix = self.prefix
        return {
            'preprocessed': prefix / 'preprocessed.interval_list',
            'annotated': prefix / 'annotated_intervals.tsv',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(multicohort)
        jobs = gcnv.prepare_intervals(self.get_job_attrs(multicohort), outputs)
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(required_stages=PrepareIntervals)
class CollectReadCounts(SequencingGroupStage):
    """
    Per-sample stage that runs CollectReadCounts to produce .counts.tsv.gz files.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'counts': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz',
            'index': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz.tbi',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        if seqgroup.cram is None:
            raise ValueError(f'No CRAM file found for {seqgroup}')

        jobs = gcnv.collect_read_counts(
            intervals_path=inputs.as_path(get_multicohort(), PrepareIntervals, 'preprocessed'),
            cram_path=seqgroup.cram,
            job_attrs=self.get_job_attrs(seqgroup),
            output_base_path=seqgroup.dataset.prefix() / 'gcnv' / seqgroup.id,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=[SetSGIDOrdering, PrepareIntervals, CollectReadCounts])
class DeterminePloidy(CohortStage):
    """
    The non-sharded cohort-wide gCNV steps after read counts have been collected:
    FilterIntervals and DetermineGermlineContigPloidy. These outputs represent
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'filtered': cohort_prefix / 'filtered.interval_list',
            'calls': cohort_prefix / 'ploidy-calls.tar.gz',
            'model': cohort_prefix / 'ploidy-model.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        prep_intervals = inputs.as_dict(get_multicohort(), PrepareIntervals)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(cohort, SetSGIDOrdering, 'sgid_order').open())
        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')
        # order those WRT the set ordering
        ordered_read_counts = [random_read_counts[seqgroup] for seqgroup in sgid_ordering]

        jobs = gcnv.filter_and_determine_ploidy(
            ploidy_priors_path=reference_path('gatk_sv/contig_ploidy_priors'),
            preprocessed_intervals_path=prep_intervals['preprocessed'],
            annotated_intervals_path=prep_intervals['annotated'],
            counts_paths=ordered_read_counts,
            job_attrs=self.get_job_attrs(cohort),
            output_paths=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=DeterminePloidy)
class UpgradePedWithInferred(CohortStage):
    """
    Don't trust the metamist pedigrees, update with inferred sexes
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'aneuploidy_samples': cohort_prefix / 'aneuploidies.txt',
            'pedigree': cohort_prefix / 'inferred_sex_pedigree.ped',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(cohort)
        ploidy_inputs = get_batch().read_input(str(inputs.as_dict(cohort, DeterminePloidy)['calls']))
        tmp_ped_path = get_batch().read_input(
            str(
                cohort.write_ped_file(self.get_stage_cohort_prefix(cohort, category='tmp') / 'pedigree.ped'),
            ),
        )
        job = gcnv.upgrade_ped_file(
            local_ped=tmp_ped_path,
            new_output=str(outputs['pedigree']),
            aneuploidies=str(outputs['aneuploidy_samples']),
            ploidy_tar=ploidy_inputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=job)  # type: ignore


@stage(required_stages=[SetSGIDOrdering, PrepareIntervals, CollectReadCounts, DeterminePloidy])
class GermlineCNV(CohortStage):
    """
    The cohort-wide GermlineCNVCaller step, sharded across genome regions.
    This is separate from the DeterminePloidy stage so that the GermlineCNVCalls
    stage can pick out this stage's sharded inputs easily.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {name: self.get_stage_cohort_prefix(cohort) / f'{name}.tar.gz' for name in gcnv.shard_basenames()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)
        determine_ploidy = inputs.as_dict(cohort, DeterminePloidy)
        prep_intervals = inputs.as_dict(get_multicohort(), PrepareIntervals)
        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(cohort, SetSGIDOrdering, 'sgid_order').open())
        # order per-sgid files WRT the set ordering
        ordered_read_counts = [random_read_counts[seqgroup] for seqgroup in sgid_ordering]

        jobs = gcnv.shard_gcnv(
            annotated_intervals_path=prep_intervals['annotated'],
            filtered_intervals_path=determine_ploidy['filtered'],
            ploidy_calls_path=determine_ploidy['calls'],
            counts_paths=ordered_read_counts,
            job_attrs=self.get_job_attrs(cohort),
            output_paths=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[SetSGIDOrdering, DeterminePloidy, GermlineCNV])
class GermlineCNVCalls(SequencingGroupStage):
    """
    Produces final individual VCF results by running PostprocessGermlineCNVCalls.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        """
        output paths here are per-SGID, but stored in the directory structure indicating the whole MCohort
        """

        # identify the cohort that contains this SGID
        this_cohort: Cohort = get_cohort_for_sgid(seqgroup.id)

        # this job runs per sample, on results with a cohort context
        # so we need to write the outputs to a cohort-specific location
        cohort_prefix = self.get_stage_cohort_prefix(this_cohort)
        return {
            'intervals': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'intervals_index': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'segments': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz',
            'segments_index': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'ratios': cohort_prefix / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(seqgroup)

        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(seqgroup.id)

        determine_ploidy = inputs.as_dict(this_cohort, DeterminePloidy)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(this_cohort, SetSGIDOrdering, 'sgid_order').open())

        jobs = gcnv.postprocess_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(this_cohort, GermlineCNV),
            sample_index=sgid_ordering.index(seqgroup.id),
            job_attrs=self.get_job_attrs(seqgroup),
            output_prefix=str(self.get_stage_cohort_prefix(this_cohort) / seqgroup.id),
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=[GermlineCNVCalls, UpgradePedWithInferred])
class TrimOffSexChromosomes(CohortStage):
    """
    Trim off sex chromosomes for gCNV VCFs where the SGID is detected to be Aneuploid
    The dependency chain here is a number of CohortStages, i.e. the MultiCohort as a whole
    isn't relevant to determining aneuploidy. As a result we're happy writing this to a Cohort-specific path
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path | str]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)

        # returning an empty dictionary might cause the pipeline setup to break?
        return_dict: dict[str, Path | str] = {
            'placeholder': str(cohort_prefix / 'placeholder.txt'),
        }

        # load up the file of aneuploidies - I don't think the pipeline supports passing an input directly here
        # so... I'm making a similar path and manually string-replacing it
        aneuploidy_file = str(cohort_prefix / 'aneuploidies.txt').replace(
            self.name,
            'UpgradePedWithInferred',
        )

        # optionally pick up aneuploid samples from the config
        aneuploid_samples: list[str] = config_retrieve(['gCNV', 'aneuploid_samples'], [])

        if (aneuploidy_path := to_path(aneuploidy_file)).exists():

            # read the identified aneuploidy samples file
            with aneuploidy_path.open() as handle:

                # iterate over the lines
                for line in handle:

                    # find the SGID
                    sgid = line.strip()

                    # could be an empty newline
                    if not sgid:
                        continue

                    aneuploid_samples.append(sgid)

        # log an expected output
        for sgid in set(aneuploid_samples):
            return_dict[sgid] = cohort_prefix / f'{sgid}.segments.vcf.bgz'

        return return_dict

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        For each of the SGIDs which are identified as aneuploid, create a version
        with the X & Y chromosomes trimmed off
        Plan to generate every file, so that the stage can be forced to re-run if needed
        """
        expected = self.expected_outputs(cohort)
        germline_calls = inputs.as_dict_by_target(GermlineCNVCalls)
        jobs = []
        sg_ids_in_cohort = cohort.get_sequencing_group_ids()
        for sgid, new_vcf in expected.items():
            if sgid == 'placeholder' or sgid not in sg_ids_in_cohort:
                continue
            sg_vcf = germline_calls[sgid]['segments']
            jobs.append(
                gcnv.trim_sex_chromosomes(
                    sgid,
                    str(sg_vcf),
                    str(new_vcf),
                    self.get_job_attrs(cohort),
                ),
            )
        return self.make_outputs(cohort, data=expected, jobs=jobs)  # type: ignore


@stage(
    required_stages=[
        TrimOffSexChromosomes,
        SetSGIDOrdering,
        GermlineCNVCalls,
        PrepareIntervals,
        UpgradePedWithInferred,
    ],
)
class GCNVJointSegmentation(CohortStage):
    """
    various config elements scavenged from https://github.com/broadinstitute/gatk/blob/cfd4d87ec29ac45a68f13a37f30101f326546b7d/scripts/cnv_cromwell_tests/germline/cnv_germline_case_scattered_workflow.json#L26
    continuing adaptation of https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/joint_call_exome_cnvs.wdl
    takes the individual VCFs and runs the joint segmentation step
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'clustered_vcf': cohort_prefix / 'JointClusteredSegments.vcf.gz',
            'clustered_vcf_idx': cohort_prefix / 'JointClusteredSegments.vcf.gz.tbi',
            'pedigree': cohort_prefix / 'pedigree.ped',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        So, this is a tricksy lil monster -
        Conducts a semi-heirarchical merge of the individual VCFs
        - First merge the segment files in blocks, to produce intermediate merges
        - Then merge those intermediate merges to produce the final result
        """

        # get the individual Segment VCFs
        cnv_vcfs = inputs.as_dict_by_target(GermlineCNVCalls)

        # and the dict of trimmed VCFs (can be empty)
        trimmed_vcfs = inputs.as_dict(cohort, TrimOffSexChromosomes)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(cohort, SetSGIDOrdering, 'sgid_order').open())

        # for each SGID, either get the sex chrom-trimmed one, or the default
        all_vcfs: list[str] = []
        for sgid in sgid_ordering:
            if sgid in trimmed_vcfs:
                get_logger().info(f'Using XY-trimmed VCF for {sgid}')
                all_vcfs.append(str(trimmed_vcfs[sgid]))
            elif sgid in cnv_vcfs:
                get_logger().warning(f'Using standard VCF for {sgid}')
                all_vcfs.append(str(cnv_vcfs[sgid]['segments']))
            else:
                raise ValueError(f'No VCF found for {sgid}')

        # get the intervals
        intervals = inputs.as_path(get_multicohort(), PrepareIntervals, 'preprocessed')

        expected_out = self.expected_outputs(cohort)

        pedigree = inputs.as_dict(cohort, UpgradePedWithInferred)['pedigree']

        jobs = gcnv.run_joint_segmentation(
            segment_vcfs=all_vcfs,
            pedigree=str(pedigree),
            intervals=str(intervals),
            tmp_prefix=str(self.get_stage_cohort_prefix(cohort, category='tmp') / 'intermediate_jointseg'),
            output_path=expected_out['clustered_vcf'],  # type: ignore
            job_attrs=self.get_job_attrs(cohort),
        )
        return self.make_outputs(cohort, data=expected_out, jobs=jobs)


@stage(
    required_stages=[SetSGIDOrdering, GCNVJointSegmentation, GermlineCNV, GermlineCNVCalls, DeterminePloidy],
)
class RecalculateClusteredQuality(SequencingGroupStage):
    """
    following joint segmentation, we need to post-process the clustered breakpoints
    this recalculates each sample's quality scores based on new breakpoints, and
    filters low QS or high AF calls
    https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/joint_call_exome_cnvs.wdl#L113

    This is done as another pass through PostprocessGermlineCNVCalls, with prior/clustered results
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:

        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(seqgroup.id)

        cohort_prefix = self.get_stage_cohort_prefix(this_cohort)

        # this job runs per sample, on results with a cohort context
        # so we need to write the outputs to a cohort-specific location
        return {
            'genotyped_intervals_vcf': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'genotyped_intervals_vcf_index': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'genotyped_segments_vcf': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz',
            'genotyped_segments_vcf_index': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'denoised_copy_ratios': cohort_prefix / f'{seqgroup.id}.ratios.tsv',
            'qc_status_file': cohort_prefix / f'{seqgroup.id}.qc_status.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(sequencing_group)

        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(sequencing_group.id)

        # get the clustered VCF from the previous stage
        joint_seg = inputs.as_dict(this_cohort, GCNVJointSegmentation)

        determine_ploidy = inputs.as_dict(this_cohort, DeterminePloidy)
        gcnv_call_inputs = inputs.as_dict(sequencing_group, GermlineCNVCalls)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(this_cohort, SetSGIDOrdering, 'sgid_order').open())

        jobs = gcnv.postprocess_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(this_cohort, GermlineCNV),
            sample_index=sgid_ordering.index(sequencing_group.id),
            job_attrs=self.get_job_attrs(sequencing_group),
            output_prefix=str(self.get_stage_cohort_prefix(this_cohort) / sequencing_group.id),
            clustered_vcf=str(joint_seg['clustered_vcf']),
            intervals_vcf=str(gcnv_call_inputs['intervals']),
            qc_file=str(expected_out['qc_status_file']),
        )
        return self.make_outputs(sequencing_group, data=expected_out, jobs=jobs)


@stage(required_stages=RecalculateClusteredQuality)
class FastCombineGCNVs(CohortStage):
    """
    Produces final multi-sample VCF results by running a merge
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        This is now explicitly continuing from multicohort work, so the output path must include
        pointers to both the MultiCohort and the Cohort
        """
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'combined_calls': cohort_prefix / 'gcnv_joint_call.vcf.bgz',
            'combined_calls_index': cohort_prefix / 'gcnv_joint_call.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        # do a slapdash bcftools merge on all input files...
        gcnv_vcfs = inputs.as_dict_by_target(RecalculateClusteredQuality)
        all_vcfs = [str(gcnv_vcfs[sgid]['genotyped_segments_vcf']) for sgid in cohort.get_sequencing_group_ids()]

        pipeline_image = get_images(['sv_pipeline_docker'])['sv_pipeline_docker']

        job_or_none = gcnv.merge_calls(
            sg_vcfs=all_vcfs,
            docker_image=pipeline_image,
            job_attrs=self.get_job_attrs(cohort),
            output_path=outputs['combined_calls'],
        )
        return self.make_outputs(cohort, data=outputs, jobs=job_or_none)


@stage(required_stages=FastCombineGCNVs, analysis_type='cnv', analysis_keys=['merged_vcf'])
class MergeCohortsgCNV(MultiCohortStage):
    """
    Takes all the per-Cohort results and merges them into a pseudocallset
    We could use Jasmine for a better merge
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        prefix = self.prefix
        return {
            'merged_vcf': prefix / 'multi_cohort_gcnv.vcf.bgz',
            'merged_vcf_index': prefix / 'multi_cohort_gcnv.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """

        Args:
            multicohort ():
            inputs (StageInput): link to FastCombineGCNVs outputs

        Returns:
            the bcftools merge job to join the Cohort-VCFs into a MultiCohort VCF
        """
        outputs = self.expected_outputs(multicohort)
        cohort_merges = inputs.as_dict_by_target(FastCombineGCNVs)
        cohort_vcfs = [str(cohort_merges[cohort.id]['combined_calls']) for cohort in multicohort.get_cohorts()]

        pipeline_image = get_images(['sv_pipeline_docker'])['sv_pipeline_docker']

        job_or_none = gcnv.merge_calls(
            sg_vcfs=cohort_vcfs,
            docker_image=pipeline_image,
            job_attrs=self.get_job_attrs(multicohort),
            output_path=outputs['merged_vcf'],
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job_or_none)


@stage(required_stages=MergeCohortsgCNV, analysis_type='cnv', analysis_keys=['annotated_vcf'])
class AnnotateCNV(MultiCohortStage):
    """
    Smaller, direct annotation using SvAnnotate
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    This is a full clone of the GATK-SV pipeline Cromwell stage, but use on a slightly
    different output. Trying to work out the best way to handle this through inheritance

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        prefix = self.prefix
        return {
            'annotated_vcf': prefix / 'merged_gcnv_annotated.vcf.bgz',
            'annotated_vcf_index': prefix / 'merged_gcnv_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        expected_out = self.expected_outputs(multicohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}
        job_or_none = queue_annotate_sv_jobs(
            multicohort=multicohort,
            prefix=self.prefix,
            input_vcf=inputs.as_dict(multicohort, MergeCohortsgCNV)['merged_vcf'],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(multicohort, data=expected_out, jobs=job_or_none)


@stage(required_stages=AnnotateCNV, analysis_type='cnv', analysis_keys=['strvctvre_vcf'])
class AnnotateCNVVcfWithStrvctvre(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        prefix = self.prefix
        return {
            'strvctvre_vcf': prefix / 'cnv_strvctvre_annotated.vcf.bgz',
            'strvctvre_vcf_index': prefix / 'cnv_strvctvre_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        strv_job = get_batch().new_job('StrVCTVRE', self.get_job_attrs() | {'tool': 'strvctvre'})

        strv_job.image(image_path('strvctvre'))
        strv_job.cpu(config_retrieve(['strvctvre_resources', 'cpu'], 2))
        strv_job.memory(config_retrieve(['strvctvre_resources', 'memory'], '20Gi'))
        strv_job.storage(config_retrieve(['strvctvre_resources', 'storage'], '20Gi'))

        strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']
        assert isinstance(strvctvre_phylop, str)
        phylop_in_batch = get_batch().read_input(strvctvre_phylop)

        input_dict = inputs.as_dict(multicohort, AnnotateCNV)
        expected_d = self.expected_outputs(multicohort)

        # read vcf and index into the batch
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']

        strv_job.declare_resource_group(output_vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # run strvctvre
        strv_job.command(f'python StrVCTVRE.py -i {input_vcf} -o temp.vcf -f vcf -p {phylop_in_batch}')
        strv_job.command(f'bgzip temp.vcf -c > {strv_job.output_vcf["vcf.bgz"]}')  # type: ignore
        strv_job.command(f'tabix {strv_job.output_vcf["vcf.bgz"]}')  # type: ignore

        get_batch().write_output(
            strv_job.output_vcf,
            str(expected_d['strvctvre_vcf']).replace('.vcf.bgz', ''),
        )
        return self.make_outputs(multicohort, data=expected_d, jobs=strv_job)


@stage(required_stages=AnnotateCNVVcfWithStrvctvre, analysis_type='cnv', analysis_keys=['mt'])
class AnnotateCohortgCNV(MultiCohortStage):
    """
    Rearrange the annotations across the cohort to suit Seqr
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {'mt': self.prefix / 'gcnv_annotated_cohort.mt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Fire up the job to ingest the cohort VCF as a MT, and rearrange the annotations
        """

        vcf_path = str(inputs.as_path(target=multicohort, stage=AnnotateCNVVcfWithStrvctvre, key='strvctvre_vcf'))
        outputs = self.expected_outputs(multicohort)

        j = get_batch().new_job('annotate gCNV cohort', self.get_job_attrs(multicohort))
        j.image(config_retrieve(['workflow', 'driver_image']))
        gencode_gz = config_retrieve(['workflow', 'gencode_gtf_file'])
        gencode_gtf_local = get_batch().read_input(str(gencode_gz))
        j.command(
            'seqr_loader_cnv '
            f'--mt_out {str(outputs["mt"])} '
            f'--checkpoint {str(self.tmp_prefix / "checkpoints")} '
            f'cohort '  # use the annotate_COHORT functionality
            f'--vcf {vcf_path} '
            f'--gencode {gencode_gtf_local} ',
        )

        return self.make_outputs(multicohort, data=outputs, jobs=j)


@stage(required_stages=AnnotateCohortgCNV, analysis_type='cnv', analysis_keys=['mt'])
class AnnotateDatasetCNV(DatasetStage):
    """
    Subset the MT to be this Dataset only, then work up all the genotype values
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to generate a matrix table
        """
        return {'mt': dataset.prefix() / 'mt' / f'gCNV-{get_workflow().output_version}-{dataset.name}.mt'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Subsets the multicohort MT to this dataset only, then brings genotype data into row annotations
        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs (StageInput): results of AnnotateCohortgCNV for this MultiCohort
        """

        mt_in = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortgCNV, key='mt')
        outputs = self.expected_outputs(dataset)

        jobs = gcnv.annotate_dataset_jobs_cnv(
            mt_path=mt_in,
            sgids=dataset.get_sequencing_group_ids(),
            out_mt_path=str(outputs['mt']),
            tmp_prefix=self.tmp_prefix / dataset.name / 'checkpoints',
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage(
    required_stages=[AnnotateDatasetCNV],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    # https://github.com/populationgenomics/metamist/issues/539
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'CNV'},
)
class MtToEsCNV(DatasetStage):
    """
    Create a Seqr index
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        index_name = f'{dataset.name}-exome-gCNV-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Freshly liberated from the clutches of DataProc
        Uses the script cpg_workflows/dataproc_scripts/mt_to_es_free_of_dataproc.py
        The script was registered in setup.py with a console entrypoint
        This requires a code version >= 1.25.14 in the worker job image to operate

        gCNV indexes have never been spotted over a GB in the wild, so we use minimum storage

        Args:
            dataset (Dataset):
            inputs ():
        """

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
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetCNV, key='mt'))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        outputs = self.expected_outputs(dataset)

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

        return self.make_outputs(dataset, data=outputs, jobs=job)
