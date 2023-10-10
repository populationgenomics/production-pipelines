"""
The second multi-sample workflow, containing the following stages:
MakeCohortVCF and AnnotateVCF
"""


from typing import Any

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_workflows.jobs.seqr_loader_sv import (
    annotate_cohort_jobs_sv,
    annotate_dataset_jobs_sv,
)
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
    make_combined_ped,
    _sv_filtered_meta,
    SV_CALLERS,
)
from cpg_workflows.workflow import (
    get_workflow,
    stage,
    Cohort,
    CohortStage,
    Dataset,
    DatasetStage,
    StageOutput,
    StageInput,
)
from cpg_workflows.jobs import ploidy_table_from_ped


@stage
class MakeCohortVcf(CohortStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    This Stage is the first in the gatk_sv_multisample_2 workflow. Using Analysis-Runner,
    include two additional config files:

    - configs/gatk_sv/use_for_all_workflows.toml
        - contains all required images and references
    - configs/gatk_sv/all_batch_names.toml
        - add all sub-cohort hashes to the batch_names list
    - A custom config with all combined SGs in workflow.only_sgs

    Once this becomes a more frequently run process, we'll have to investigate whether
    we want to include ALL prior batches, or just the batches which were run in the
    earlier pipeline steps this time around. That may depend on how much these stages
    cost, and whether generating the final cross-batch VCF/ES-index is possible from
    several component joint calls, or we require one single joint call across all SGs
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """create output paths"""
        return {
            'vcf': self.prefix / 'cleaned.vcf.gz',
            'vcf_index': self.prefix / 'cleaned.vcf.gz.tbi',
            'vcf_qc': self.prefix / 'cleaned_SV_VCF_QC_output.tar.gz',
            'metrics_file_makecohortvcf': self.prefix / 'metrics.tsv',
            # if merge_intermediate_vcfs is enabled
            # 'cluster_vcf': '.combine_batches.vcf.gz',
            # 'cluster_vcf_index': '.combine_batches.vcf.gz.tbi',
            # 'complex_resolve_vcf': '.complex_resolve.vcf.gz',
            # 'complex_resolve_vcf_index': '.complex_resolve.vcf.gz.tbi',
            # 'complex_genotype_vcf': '.complex_genotype.vcf.gz',
            # 'complex_genotype_vcf_index': '.complex_genotype.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        This is a little bit spicy. Instead of taking a direct dependency on the
        previous stage, we use the output hash to find all the previous batches

        Replacing this nasty mess with nested/fancy cohorts would be ace
        """

        batch_names = get_config()['workflow']['batch_names']
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        pesr_vcfs = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'pesr.vcf.gz'
            for batch_name in batch_names
        ]
        depth_vcfs = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'depth.vcf.gz'
            for batch_name in batch_names
        ]
        sr_pass = [
            batch_prefix
            / batch_name
            / 'GenotypeBatch'
            / 'genotype_SR_part2_bothside_pass.txt'
            for batch_name in batch_names
        ]
        sr_fail = [
            batch_prefix
            / batch_name
            / 'GenotypeBatch'
            / 'genotype_SR_part2_background_fail.txt'
            for batch_name in batch_names
        ]
        depth_depth_cutoff = [
            batch_prefix / batch_name / 'GenotypeBatch' / 'depth.depth_sepcutoff.txt'
            for batch_name in batch_names
        ]
        filter_batch_cutoffs = [
            batch_prefix / batch_name / 'FilterBatch' / 'cutoffs'
            for batch_name in batch_names
        ]
        bincov_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'GatherBatchEvidence.RD.txt.gz'
            for batch_name in batch_names
        ]
        disc_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'GatherBatchEvidence.pe.txt.gz'
            for batch_name in batch_names
        ]
        median_cov_files = [
            batch_prefix
            / batch_name
            / 'GatherBatchEvidence'
            / 'medianCov.transposed.bed'
            for batch_name in batch_names
        ]

        input_dict: dict[str, Any] = {
            'cohort_name': cohort.name,
            'batches': batch_names,
            'ped_file': make_combined_ped(cohort, self.prefix),
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
            ]
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
            ]
        )
        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=MakeCohortVcf)
class FormatVcfForGatk(CohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'gatk_formatted_vcf': self.prefix / 'gatk_formatted.vcf.gz',
            'gatk_formatted_vcf_index': self.prefix / 'gatk_formatted.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Args:
            cohort (Cohort): cohort of all samples (across several sub-cohort batches)
            inputs (StageInput): access to prior inputs
        """

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'vcf': inputs.as_dict(cohort, MakeCohortVcf)['vcf'],
            'ped_file': make_combined_ped(cohort, self.prefix),
        }
        input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=MakeCohortVcf)
class JoinRawCalls(CohortStage):
    """
    Joins all individually clustered caller results
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'joined_raw_calls_vcf': self.prefix / 'raw_clustered_calls.vcf.gz',
            'joined_raw_calls_vcf_index': self.prefix
            / 'raw_clustered_calls.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """ """

        input_dict: dict[str, Any] = {
            'prefix': cohort.name,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'reference_fasta': get_fasta(),
            'reference_fasta_fai': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
        }
        input_dict |= get_images(
            ['gatk_docker', 'sv_pipeline_docker', 'sv_base_mini_docker']
        )
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        # add all clustered _caller_ files, plus indices
        batch_names = get_config()['workflow']['batch_names']
        batch_prefix = cohort.analysis_dataset.prefix() / 'gatk_sv'
        for caller in SV_CALLERS + ['depth']:
            input_dict[f'clustered_{caller}_vcfs'] = [
                batch_prefix
                / batch_name
                / 'ClusterBatch'
                / f'clustered-{caller}.vcf.gz'
                for batch_name in batch_names
            ]
            input_dict[f'clustered_{caller}_vcf_indexes'] = [
                batch_prefix
                / batch_name
                / 'ClusterBatch'
                / f'clustered-{caller}.vcf.gz.tbi'
                for batch_name in batch_names
            ]
        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=[JoinRawCalls, FormatVcfForGatk])
class SVConcordance(CohortStage):
    """
    Takes the clean VCF and reformat for GATK intake
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'concordance_vcf': self.prefix / 'sv_concordance.vcf.gz',
            'concordance_vcf_index': self.prefix / 'sv_concordance.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV concordance
        """
        raw_calls = inputs.as_dict(cohort, JoinRawCalls)
        format_vcf = inputs.as_dict(cohort, FormatVcfForGatk)

        input_dict: dict[str, Any] = {
            'output_prefix': cohort.name,
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'eval_vcf': format_vcf['gatk_formatted_vcf'],
            'truth_vcf': raw_calls['joined_raw_calls_vcf'],
        }
        input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
        input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=SVConcordance)
class GeneratePloidyTable(CohortStage):
    """
    Quick PythonJob to generate a ploidy table
    Calls a homebrewed version of this table generator:
    github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/ploidy_table_from_ped.py
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        only one output, the ploidy table
        """

        return {'ploidy_table': self.prefix / 'ploidy_table.txt'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        py_job = get_batch().new_python_job('create_ploidy_table')
        py_job.image(get_config()['workflow']['driver_image'])

        ped_path = make_combined_ped(cohort, self.prefix)
        contig_path = get_references(['primary_contigs_list'])['primary_contigs_list']

        expected_d = self.expected_outputs(cohort)
        py_job.call(
            ploidy_table_from_ped.generate_ploidy_table,
            ped_path,
            contig_path,
            expected_d['ploidy_table'],
        )

        return self.make_outputs(cohort, data=expected_d, jobs=py_job)


@stage(
    required_stages=[GeneratePloidyTable, SVConcordance],
    analysis_type='sv',
    analysis_keys=['filtered_vcf'],
    update_analysis_meta=_sv_filtered_meta,
)
class FilterGenotypes(CohortStage):
    """
    Steps required to post-filter called genotypes
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        create dictionary of names -> output paths
        """

        return {
            'filtered_vcf': self.prefix / 'filtered.vcf.gz',
            'filtered_vcf_index': self.prefix / 'filtered.vcf.gz.tbi',
            'main_vcf_qc_tarball': self.prefix / 'filtered_SV_VCF_QC_output.tar.gz',
            'unfiltered_recalibrated_vcf': self.prefix
            / 'unfiltered_recalibrated.vcf.gz',
            'unfiltered_recalibrated_vcf_index': self.prefix
            / 'unfiltered_recalibrated.vcf.gz.tbi',
            # 'vcf_optimization_table': self.prefix / 'vcf_optimization_table.tsv.gz',
            # 'sl_cutoff_qc_tarball': self.prefix / 'sl_cutoff_SV_VCF_QC_output.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        input_dict = {
            'output_prefix': cohort.name,
            'vcf': inputs.as_dict(cohort, SVConcordance)['concordance_vcf'],
            'ploidy_table': inputs.as_dict(cohort, GeneratePloidyTable)['ploidy_table'],
            'ped_file': make_combined_ped(cohort, self.prefix),
            'fmax_beta': get_config()['references']['gatk_sv'].get('fmax_beta', 0.3),
            'recalibrate_gq_args': get_config()['references']['gatk_sv'].get(
                'recalibrate_gq_args'
            ),
            'sl_filter_args': get_config()['references']['gatk_sv'].get(
                'sl_filter_args'
            ),
        }
        assert input_dict['recalibrate_gq_args']
        assert input_dict['sl_filter_args']

        input_dict |= get_images(
            ['linux_docker', 'sv_base_mini_docker', 'sv_pipeline_docker']
        )
        # use a non-standard GATK image containing required filtering tool
        input_dict['gatk_docker'] = get_images(['gq_recalibrator_docker'])[
            'gq_recalibrator_docker'
        ]
        input_dict |= get_references(
            [
                {'gq_recalibrator_model_file': 'aou_filtering_model'},
                'primary_contigs_fai',
            ]
        )

        # something a little trickier - we need to get various genome tracks
        input_dict['genome_tracks'] = list(
            get_references(
                get_config()['references']['gatk_sv'].get('genome_tracks', [])
            ).values()
        )

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


def _sv_annotated_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Callable, add meta[type] to custom analysis object
    """
    return {'type': 'gatk-sv-annotated'}


@stage(
    required_stages=FilterGenotypes,
    analysis_type='sv',
    analysis_keys=['annotated_vcf'],
    update_analysis_meta=_sv_annotated_meta,
)
class AnnotateVcf(CohortStage):
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
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        return {
            'annotated_vcf': self.prefix / 'filtered_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'filtered_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        input_dict: dict[str, Any] = {
            'vcf': inputs.as_dict(cohort, FilterGenotypes)['filtered_vcf'],
            'prefix': cohort.name,
            'ped_file': make_combined_ped(cohort, self.prefix),
            'sv_per_shard': 5000,
            'population': get_config()['references']['gatk_sv'].get(
                'external_af_population'
            ),
            'ref_prefix': get_config()['references']['gatk_sv'].get(
                'external_af_ref_bed_prefix'
            ),
            'use_hail': False,
        }

        input_dict |= get_references(
            [
                'noncoding_bed',
                'protein_coding_gtf',
                {'ref_bed': 'external_af_ref_bed'},
                {'contig_list': 'primary_contigs_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            ['sv_pipeline_docker', 'sv_base_mini_docker', 'gatk_docker']
        )
        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=AnnotateVcf)
class AnnotateCohortSv(CohortStage):
    """
    What do we want?! SV Data in Seqr!
    When do we want it?! Now!

    First step to transform annotated SV callset data into a seqr ready format
    Rearrange all the annotations
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Expected to write a matrix table.
        """
        return {
            'tmp_prefix': str(self.tmp_prefix),
            'mt': self.prefix / 'cohort_sv.mt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """

        vcf_path = inputs.as_path(target=cohort, stage=AnnotateVcf, key='annotated_vcf')
        checkpoint_prefix = (
            to_path(self.expected_outputs(cohort)['tmp_prefix']) / 'checkpoints'
        )

        job = annotate_cohort_jobs_sv(
            b=get_batch(),
            vcf_path=vcf_path,
            out_mt_path=self.expected_outputs(cohort)['mt'],
            checkpoint_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(cohort),
            depends_on=inputs.get_jobs(cohort),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=job)


def _update_sv_dataset_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to custom analysis object
    """
    return {'type': 'annotated-sv-dataset-callset'}


@stage(
    required_stages=AnnotateCohortSv,
    analysis_type='sv',
    analysis_keys=['mt'],
    update_analysis_meta=_update_sv_dataset_meta,
)
class AnnotateDatasetSv(DatasetStage):
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
            'mt': (
                dataset.prefix()
                / 'mt'
                / f'{get_workflow().output_version}-{dataset.name}.mt'
            ),
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Whether Dataproc or not, this Stage subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """

        assert dataset.cohort
        mt_path = inputs.as_path(
            target=dataset.cohort, stage=AnnotateCohortSv, key='mt'
        )

        checkpoint_prefix = (
            to_path(self.expected_outputs(dataset)['tmp_prefix']) / 'checkpoints'
        )

        jobs = annotate_dataset_jobs_sv(
            b=get_batch(),
            mt_path=mt_path,
            sgids=dataset.get_sequencing_group_ids(),
            out_mt_path=self.expected_outputs(dataset)['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
            depends_on=inputs.get_jobs(dataset)
        )

        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


def _gatk_sv_index_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Add meta.type to custom analysis object
    https://github.com/populationgenomics/metamist/issues/539
    """
    return {'type': 'gatk-sv-index', 'seqr-data-type': 'gatk-sv'}


@stage(
    required_stages=[AnnotateDatasetSv],
    analysis_type='es-index',  # specific type of es index
    analysis_keys=['index_name'],
    update_analysis_meta=_gatk_sv_index_meta,
)
class MtToEsSv(DatasetStage):
    """
    Create a Seqr index
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = get_config()['workflow']['sequencing_type']
        index_name = f'{dataset.name}-{sequencing_type}-SV-{get_workflow().run_timestamp}'.lower()
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

        dataset_mt_path = inputs.as_path(
            target=dataset, stage=AnnotateDatasetSv, key='mt'
        )
        index_name = self.expected_outputs(dataset)['index_name']
        done_flag_path = self.expected_outputs(dataset)['done_flag']

        if 'elasticsearch' not in get_config():
            raise ValueError(
                f'"elasticsearch" section is not defined in config, cannot create '
                f'Elasticsearch index for dataset {dataset}'
            )

        from analysis_runner import dataproc

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
                init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
            )
        j._preemptible = False
        j.attributes = (j.attributes or {}) | {'tool': 'hailctl dataproc'}
        jobs = [j]
        return self.make_outputs(dataset, data=index_name, jobs=jobs)
