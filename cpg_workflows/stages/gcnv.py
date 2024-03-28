"""
Stages that implement GATK-gCNV.
"""

import json

from cpg_utils import Path, to_path
from cpg_utils.config import AR_GUID_NAME, get_config, try_get_ar_guid
from cpg_utils.hail_batch import get_batch, image_path, query_command, reference_path
from cpg_workflows.inputs import get_cohort
from cpg_workflows.jobs import gcnv
from cpg_workflows.query_modules import seqr_loader_cnv
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    get_images,
    get_references,
    queue_annotate_sv_jobs,
)
from cpg_workflows.targets import Cohort, SequencingGroup
from cpg_workflows.workflow import (
    CohortStage,
    Dataset,
    DatasetStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)


def _gcnv_annotated_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Callable, adds custom analysis object meta attribute
    """
    return {'type': 'gCNV-annotated'}


@stage
class SetSGIDOrdering(CohortStage):
    """
    Set the order of the sequencing groups in the cohort
    Push this to a file _now_, and use it later
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {'sgid_order': self.prefix / 'sgid_order.json'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        sorted_sgids = sorted(cohort.get_sequencing_group_ids())
        with self.expected_outputs(cohort)['sgid_order'].open('w') as f_handler:
            json.dump(sorted_sgids, f_handler, indent=2)
        return self.make_outputs(cohort)


@stage
class PrepareIntervals(CohortStage):
    """
    Interval preparation steps that don't require the sample read counts:
    PreprocessIntervals and AnnotateIntervals.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'preprocessed': self.prefix / 'preprocessed.interval_list',
            'annotated': self.prefix / 'annotated_intervals.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)
        jobs = gcnv.prepare_intervals(self.get_job_attrs(cohort), outputs)
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


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
            intervals_path=inputs.as_path(get_cohort(), PrepareIntervals, 'preprocessed'),
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
    intermediate results for the cohort as a whole, so are written to tmp_prefix.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'filtered': self.prefix / 'filtered.interval_list',
            'calls': self.prefix / 'ploidy-calls.tar.gz',
            'model': self.prefix / 'ploidy-model.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        prep_intervals = inputs.as_dict(cohort, PrepareIntervals)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(cohort, SetSGIDOrdering, 'sgid_order').open())
        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')
        # order those WRT the set ordering
        ordered_read_counts = [random_read_counts[seqgroup] for seqgroup in sgid_ordering]

        jobs = gcnv.filter_and_determine_ploidy(
            ploidy_priors_path=str(reference_path('gatk_sv/contig_ploidy_priors')),
            preprocessed_intervals_path=prep_intervals['preprocessed'],
            annotated_intervals_path=prep_intervals['annotated'],
            counts_paths=ordered_read_counts,
            job_attrs=self.get_job_attrs(cohort),
            output_paths=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[SetSGIDOrdering, PrepareIntervals, CollectReadCounts, DeterminePloidy])
class GermlineCNV(CohortStage):
    """
    The cohort-wide GermlineCNVCaller step, sharded across genome regions.
    This is separate from the DeterminePloidy stage so that the GermlineCNVCalls
    stage can pick out this stage's sharded inputs easily.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {name: self.prefix / f'{name}.tar.gz' for name in gcnv.shard_basenames()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)
        determine_ploidy = inputs.as_dict(cohort, DeterminePloidy)
        prep_intervals = inputs.as_dict(cohort, PrepareIntervals)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(cohort, SetSGIDOrdering, 'sgid_order').open())
        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')
        # order those WRT the set ordering
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
        return {
            'intervals': self.prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'intervals_index': self.prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'segments': self.prefix / f'{seqgroup.id}.segments.vcf.gz',
            'segments_index': self.prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'ratios': self.prefix / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(seqgroup)
        determine_ploidy = inputs.as_dict(get_cohort(), DeterminePloidy)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(get_cohort(), SetSGIDOrdering, 'sgid_order').open())

        jobs = gcnv.postprocess_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(get_cohort(), GermlineCNV),
            sample_index=sgid_ordering.index(seqgroup.id),
            job_attrs=self.get_job_attrs(seqgroup),
            output_prefix=str(self.prefix / seqgroup.id),
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=[SetSGIDOrdering, GermlineCNVCalls, PrepareIntervals])
class GCNVJointSegmentation(CohortStage):
    """
    various config elements scavenged from https://github.com/broadinstitute/gatk/blob/cfd4d87ec29ac45a68f13a37f30101f326546b7d/scripts/cnv_cromwell_tests/germline/cnv_germline_case_scattered_workflow.json#L26
    continuing adaptation of https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/joint_call_exome_cnvs.wdl
    takes the individual VCFs and runs the joint segmentation step
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, str | Path]:
        return {
            'clustered_vcf': self.prefix / 'JointClusteredSegments.vcf.gz',
            'clustered_vcf_idx': self.prefix / 'JointClusteredSegments.vcf.gz.tbi',
            'pedigree': self.tmp_prefix / 'pedigree.ped',
            'tmp_prefix': str(self.tmp_prefix / 'intermediate_jointseg'),
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        So, this is a tricksy lil monster -
        Conducts a semi-heirarchical merge of the individual VCFs
        - First merge the segment files in blocks, to produce intermediate merges
        - Then merge those intermediate merges to produce the final result
        """

        # get the list of individual Segment VCFs
        cnv_vcfs = inputs.as_dict_by_target(GermlineCNVCalls)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(get_cohort(), SetSGIDOrdering, 'sgid_order').open())

        all_vcfs = [str(cnv_vcfs[sgid]['segments']) for sgid in sgid_ordering]

        # get the intervals
        intervals = inputs.as_path(cohort, PrepareIntervals, 'preprocessed')

        expected_out = self.expected_outputs(cohort)

        ped_path = cohort.write_ped_file(expected_out['pedigree'])

        jobs = gcnv.run_joint_segmentation(
            segment_vcfs=all_vcfs,
            pedigree=str(ped_path),
            intervals=str(intervals),
            tmp_prefix=expected_out['tmp_prefix'],
            output_path=expected_out['clustered_vcf'],
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

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        return {
            'genotyped_intervals_vcf': self.prefix / f'{sequencing_group.id}.intervals.vcf.gz',
            'genotyped_intervals_vcf_index': self.prefix / f'{sequencing_group.id}.intervals.vcf.gz.tbi',
            'genotyped_segments_vcf': self.prefix / f'{sequencing_group.id}.segments.vcf.gz',
            'genotyped_segments_vcf_index': self.prefix / f'{sequencing_group.id}.segments.vcf.gz.tbi',
            'denoised_copy_ratios': self.prefix / f'{sequencing_group.id}.ratios.tsv',
            'qc_status_file': self.prefix / f'{sequencing_group.id}.qc_status.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(sequencing_group)

        # get the clustered VCF from the previous stage
        joint_seg = inputs.as_dict(get_cohort(), GCNVJointSegmentation)

        determine_ploidy = inputs.as_dict(get_cohort(), DeterminePloidy)
        gcnv_call_inputs = inputs.as_dict(sequencing_group, GermlineCNVCalls)

        # pull the json file with the sgid ordering
        sgid_ordering = json.load(inputs.as_path(get_cohort(), SetSGIDOrdering, 'sgid_order').open())

        jobs = gcnv.postprocess_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(get_cohort(), GermlineCNV),
            sample_index=sgid_ordering.index(sequencing_group.id),
            job_attrs=self.get_job_attrs(sequencing_group),
            output_prefix=str(self.prefix / sequencing_group.id),
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
        return {
            'combined_calls': self.prefix / 'gcnv_joint_call.vcf.bgz',
            'combined_calls_index': self.prefix / 'gcnv_joint_call.vcf.bgz.tbi',
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


@stage(
    required_stages=FastCombineGCNVs,
    analysis_type='cnv',
    analysis_keys=['annotated_vcf'],
    update_analysis_meta=_gcnv_annotated_meta,
)
class AnnotateCNV(CohortStage):
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

    def expected_outputs(self, cohort: Cohort) -> dict:
        return {
            'annotated_vcf': self.prefix / 'gcnv_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'gcnv_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        expected_out = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        job_or_none = queue_annotate_sv_jobs(
            cohort=cohort,
            cohort_prefix=self.prefix,
            input_vcf=inputs.as_dict(cohort, FastCombineGCNVs)['combined_calls'],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_out, jobs=job_or_none)


def _gcnv_srvctvre_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Callable, adds custom analysis object meta attribute
    """
    return {'type': 'gCNV-STRVCTCRE-annotated'}


@stage(
    required_stages=AnnotateCNV,
    analysis_type='cnv',
    analysis_keys=['strvctvre_vcf'],
    update_analysis_meta=_gcnv_srvctvre_meta,
)
class AnnotateCNVVcfWithStrvctvre(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.prefix / 'cnv_strvctvre_annotated.vcf.bgz',
            'strvctvre_vcf_index': self.prefix / 'cnv_strvctvre_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        strv_job = get_batch().new_job('StrVCTVRE', self.get_job_attrs() | {'tool': 'strvctvre'})

        strv_job.image(image_path('strvctvre'))
        strv_job.storage('20Gi')
        strv_job.memory('8Gi')

        strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']
        phylop_in_batch = get_batch().read_input(strvctvre_phylop)

        input_dict = inputs.as_dict(cohort, AnnotateCNV)
        expected_d = self.expected_outputs(cohort)

        # read vcf and index into the batch
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']

        strv_job.declare_resource_group(
            output_vcf={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # run strvctvre
        strv_job.command(f'python StrVCTVRE.py -i {input_vcf} -o temp.vcf -f vcf -p {phylop_in_batch}')
        strv_job.command(f'bgzip temp.vcf -c > {strv_job.output_vcf["vcf.bgz"]}')
        strv_job.command(f'tabix {strv_job.output_vcf["vcf.bgz"]}')

        get_batch().write_output(
            strv_job.output_vcf,
            str(expected_d['strvctvre_vcf']).replace('.vcf.bgz', ''),
        )
        return self.make_outputs(cohort, data=expected_d, jobs=strv_job)


@stage(
    required_stages=AnnotateCNVVcfWithStrvctvre,
    analysis_type='cnv',
    analysis_keys=['mt'],
    update_analysis_meta=_gcnv_srvctvre_meta,
)
class AnnotateCohortgCNV(CohortStage):
    """
    Rearrange the annotations across the cohort to suit Seqr
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path | str]:
        # convert temp path to str to avoid checking existence
        return {
            'tmp_prefix': str(self.tmp_prefix),
            'mt': self.prefix / 'gcnv_annotated_cohort.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Fire up the job to ingest the cohort VCF as a MT, and rearrange the annotations
        """

        vcf_path = inputs.as_path(target=cohort, stage=AnnotateCNVVcfWithStrvctvre, key='strvctvre_vcf')

        checkpoint_prefix = to_path(self.expected_outputs(cohort)['tmp_prefix']) / 'checkpoints'
        j = get_batch().new_job('annotate gCNV cohort', self.get_job_attrs(cohort))
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                seqr_loader_cnv,
                seqr_loader_cnv.annotate_cohort_gcnv.__name__,
                str(vcf_path),
                str(self.expected_outputs(cohort)['mt']),
                str(checkpoint_prefix),
                setup_gcp=True,
            ),
        )

        # todo is this necessary?
        if depends_on := inputs.get_jobs(cohort):
            j.depends_on(*depends_on)

        return self.make_outputs(
            cohort,
            data=self.expected_outputs(cohort),
            jobs=j,
        )


@stage(required_stages=AnnotateCohortgCNV, analysis_type='cnv', analysis_keys=['mt'])
class AnnotateDatasetCNV(DatasetStage):
    """
    Subset the MT to be this Dataset only
    Then work up all the genotype values
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Expected to generate a matrix table
        """
        return {
            'tmp_prefix': str(self.tmp_prefix / dataset.name),
            'mt': (dataset.prefix() / 'mt' / f'gCNV-{get_workflow().output_version}-{dataset.name}.mt'),
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations
        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """

        assert dataset.cohort
        mt_path = inputs.as_path(target=dataset.cohort, stage=AnnotateCohortgCNV, key='mt')

        checkpoint_prefix = to_path(self.expected_outputs(dataset)['tmp_prefix']) / 'checkpoints'

        jobs = gcnv.annotate_dataset_jobs_cnv(
            mt_path=mt_path,
            sgids=dataset.get_sequencing_group_ids(),
            out_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset),  # todo is this necessary?
        )

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)


def _gatk_gcnv_index_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Add meta.type to custom analysis object
    https://github.com/populationgenomics/metamist/issues/539
    """
    return {'seqr-dataset-type': 'CNV'}


@stage(
    required_stages=[AnnotateDatasetCNV],
    analysis_type='es-index',  # specific type of es index
    analysis_keys=['index_name'],
    update_analysis_meta=_gatk_gcnv_index_meta,
)
class MtToEsCNV(DatasetStage):
    """
    Create a Seqr index
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = get_config()['workflow']['sequencing_type']
        index_name = f'{dataset.name}-{sequencing_type}-gCNV-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses analysis-runner's dataproc helper to run a hail query script
        """
        if (
            es_datasets := get_config()['workflow'].get('create_es_index_for_datasets')
        ) and dataset.name not in es_datasets:
            # Skipping dataset that wasn't explicitly requested to upload to ES
            return self.make_outputs(dataset)

        dataset_mt_path = inputs.as_path(target=dataset, stage=AnnotateDatasetCNV, key='mt')
        index_name = self.expected_outputs(dataset)['index_name']
        done_flag_path = self.expected_outputs(dataset)['done_flag']

        if 'elasticsearch' not in get_config():
            raise ValueError(
                f'"elasticsearch" section is not defined in config, cannot create '
                f'Elasticsearch index for dataset {dataset}',
            )

        from analysis_runner import dataproc
        from cpg_workflows.stages.seqr_loader import es_password

        # transformation is the same, just use the same methods file?
        script = (
            f'cpg_workflows/dataproc_scripts/mt_to_es.py '
            f'--mt-path {dataset_mt_path} '
            f'--es-index {index_name} '
            f'--done-flag-path {done_flag_path} '
            f'--es-password {es_password()}'
        )
        pyfiles = ['seqr-loading-pipelines/hail_scripts']
        job_name = f'{dataset.name}: create ES index'

        if cluster_id := get_config()['hail'].get('dataproc', {}).get('cluster_id'):
            # noinspection PyProtectedMember
            j = dataproc._add_submit_job(
                batch=get_batch(),
                cluster_id=cluster_id,
                script=script,
                pyfiles=pyfiles,
                job_name=job_name,
                region='australia-southeast1',
            )
        else:
            j = dataproc.hail_dataproc_job(
                get_batch(),
                script,
                max_age='48h',
                packages=[
                    'cpg_workflows',
                    'elasticsearch==8.*',
                    'google',
                    'fsspec',
                    'gcloud',
                ],
                num_workers=2,
                num_secondary_workers=0,
                job_name=job_name,
                depends_on=inputs.get_jobs(dataset),
                scopes=['cloud-platform'],
                pyfiles=pyfiles,
            )
        j._preemptible = False
        j.attributes = (j.attributes or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]
        return self.make_outputs(dataset, data=index_name, jobs=jobs)
