"""
Driver for computing structural variants from GATK-SV from WGS data.
"""
from os.path import join
from typing import Any

from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, reference_path, image_path, genome_build
from cpg_workflows.batch import make_job_name, Batch, get_batch
from cpg_workflows.workflow import (
    stage,
    SampleStage,
    StageOutput,
    DatasetStage,
    StageInput,
    Sample,
    Dataset,
)

GATK_SV_COMMIT = 'b59a8b070da48ceed475814117787cf80cace170'
SV_CALLERS = ['manta', 'wham', 'scramble']


def get_images(keys: list[str]) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.
    """
    return {k: image_path(k) for k in get_config()['images'].keys() if k in keys}


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        if isinstance(key, dict):
            key, ref_d_key = list(key.items())[0]
        else:
            ref_d_key = key
        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = str(reference_path(f'gatk_sv/{ref_d_key}'))
        except KeyError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))

    return res


def add_gatk_sv_jobs(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path],
    sample_id: str | None = None,
    driver_image: str | None = None,
) -> list[Job]:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sample_id:
        output_prefix = join(output_prefix, sample_id)

    outputs_to_collect = dict()
    for key in expected_out_dict.keys():
        outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

    driver_image = driver_image or image_path('cpg_workflows')

    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset.name)
    submit_j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=get_config()['workflow']['dataset'],
        access_level=get_config()['workflow']['access_level'],
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict={
            f'{wfl_name}.{k}': str(v) if isinstance(v, Path) else v
            for k, v in input_dict.items()
        },
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(driver_image)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        cmds.append(f'gsutil cp $(cat {resource}) {out_path}')
    copy_j.command(command(cmds, setup_gcp=True))
    return [submit_j, copy_j]


def get_ref_panel(keys: list[str] | None = None) -> dict:
    return {
        k: v
        for k, v in {
            'ref_panel_samples': get_config()['sv_ref_panel']['ref_panel_samples'],
            'ref_ped_file': str(reference_path('gatk_sv/ped_file')),
            'ref_panel_vcf': str(reference_path('gatk_sv/clean_vcf')),
            'qc_definitions': str(reference_path('gatk_sv/qc_definitions')),
            'ref_panel_bincov_matrix': str(
                reference_path('gatk_sv/ref_panel_bincov_matrix')
            ),
            'contig_ploidy_model_tar': str(
                reference_path('gatk_sv/contig_ploidy_model_tar')
            ),
            'gcnv_model_tars': [
                str(reference_path('gatk_sv/model_tar_tmpl')).format(shard=i)
                for i in range(get_config()['sv_ref_panel']['model_tar_cnt'])
            ],
            'ref_panel_PE_files': [
                str(reference_path('gatk_sv/ref_panel_PE_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SR_files': [
                str(reference_path('gatk_sv/ref_panel_SR_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SD_files': [
                str(reference_path('gatk_sv/ref_panel_SD_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
        }.items()
        if not keys or k in keys
    }


@stage
class GatherSampleEvidence(SampleStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Expected to produce coverage counts, a VCF for each variant caller,
        and a txt for each type of SV evidence (SR, PE, SD).
        """
        d: dict[str, Path] = dict()

        # Coverage counts
        ending_by_key = {'coverage_counts': 'coverage_counts.tsv.gz'}

        # Caller's VCFs
        for caller in SV_CALLERS:
            ending_by_key[f'{caller}_vcf'] = f'{caller}.vcf.gz'
            ending_by_key[f'{caller}_index'] = f'{caller}.vcf.gz.tbi'

        # SV evidence
        # split reads:
        ending_by_key['sr'] = 'sr.txt.gz'
        ending_by_key['sr_index'] = 'sr.txt.gz.tbi'
        # discordant paired end reads:
        ending_by_key['pe'] = 'pe.txt.gz'
        ending_by_key['pe_index'] = 'pe.txt.gz.tbi'
        # site depth:
        ending_by_key['sd'] = 'sd.txt.gz'
        ending_by_key['sd_index'] = 'sd.txt.gz.tbi'

        for key, ending in ending_by_key.items():
            stage_name = self.name.lower()
            fname = f'{sample.id}.{ending}'  # the dot separator is critical here!!
            # it's critical to separate the ending with a dot above: `*.sr.txt.gz`,
            # `*.pe.txt.gz`, and `*.sd.txt.gz`. These files are passed to
            # `gatk PrintSVEvidence`, that would determine file format based on the
            # file name. It would strongly expect the files to end exactly with either
            # `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
            # "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
            # no suitable codecs found".
            d[key] = sample.dataset.prefix() / 'gatk_sv' / stage_name / fname
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """Add jobs to batch"""
        assert sample.cram, sample

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(sample.cram),
            'bam_or_cram_index': str(sample.cram) + '.crai',
            'sample_id': sample.id,
            # This option forces CRAM localisation, otherwise it would be passed to
            # samtools as a URL (in CramToBam.wdl) and it would fail to read it as
            # GCS_OAUTH_TOKEN is not set.
            'requester_pays_crams': True,
        }

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'samtools_cloud_docker',
                'gatk_docker',
                'genomes_in_the_cloud_docker',
                'cloud_sdk_docker',
                'wham_docker',
                'manta_docker',
                'scramble_docker',
                'sv_pipeline_base_docker',
                'gatk_docker_pesr_override',
            ]
        )
        input_dict |= get_references(
            [
                'primary_contigs_fai',
                'primary_contigs_list',
                'preprocessed_intervals',
                'manta_region_bed',
                'wham_include_list_bed_file',
                {'sd_locs_vcf': 'dbsnp_vcf'},
            ]
        )
        ref_fasta = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
        input_dict |= {
            'reference_fasta': str(ref_fasta),
            'reference_index': str(ref_fasta) + '.fai',
            'reference_dict': str(ref_fasta.with_suffix('.dict')),
        }
        input_dict['reference_version'] = '38'

        expected_d = self.expected_outputs(sample)

        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=sample.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
        )
        return self.make_outputs(sample, data=expected_d, jobs=jobs)


@stage(required_stages=GatherSampleEvidence)
class EvidenceQC(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#evidenceqc
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'ploidy_matrix': 'ploidy_matrix.bed.gz',
            'ploidy_plots': 'ploidy_plots.tar.gz',
            'WGD_dist': 'WGD_score_distributions.pdf',
            'WGD_matrix': 'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores': 'WGD_scores.txt.gz',
            'bincov_matrix': 'RD.txt.gz',
            'bincov_matrix_index': 'RD.txt.gz.tbi',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        for key, fname in fname_by_key.items():
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        d = inputs.as_dict_by_target(GatherSampleEvidence)

        sids = dataset.get_sample_ids()

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'run_vcf_qc': True,  # generates <caller>_qc_low/<caller>_qc_high
            'counts': [str(d[sid]['coverage_counts']) for sid in sids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(d[sid][f'{caller}_vcf']) for sid in sids
            ]

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
            ]
        )

        input_dict |= get_references(
            [
                'genome_file',
                'wgd_scoring_mask',
            ]
        )

        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


def make_combined_ped(dataset: Dataset) -> Path:
    """
    Create cohort + ref panel PED.
    Since the PED file must include reference panel samples as well, so concatenating
    the `dataset`'s PED file with the reference panel PED file.
    """
    combined_ped_path = (
        dataset.tmp_prefix()
        / 'gatk-sv'
        / 'ped_with_ref_panel'
        / f'{dataset.alignment_inputs_hash()}.ped'
    )
    with combined_ped_path.open('w') as out:
        with dataset.write_ped_file().open() as f:
            out.write(f.read())
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with reference_path('gatk_sv/ped_file').open() as f:
            out.write(f.read())
    return combined_ped_path


@stage(required_stages=[GatherSampleEvidence, EvidenceQC])
class GatherBatchEvidence(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#gather-batch-evidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherBatchEvidence.wdl
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        d: dict[str, Path] = dict()
        ending_by_key = {
            'cnmops_dup': '.DUP.header.bed.gz',
            'cnmops_dup_index': '.DUP.header.bed.gz.tbi',
            'cnmops_del': '.DEL.header.bed.gz',
            'cnmops_del_index': '.DEL.header.bed.gz.tbi',
            'cnmops_large_del': '.DEL.large.bed.gz',
            'cnmops_large_del_index': '.DEL.large.bed.gz.tb',
            'cnmops_large_dup': '.DUP.large.bed.gz',
            'cnmops_large_dup_index': '.DUP.large.bed.gz.tbi',
            'merged_SR': '.sr.txt.gz',
            'merged_SR_index': '.sr.txt.gz.tbi',
            'merged_PE': '.pe.txt.gz',
            'merged_PE_index': '.pe.txt.gz.tbi',
            'merged_BAF': '.baf.txt.gz',
            'merged_BAF_index': '.baf.txt.gz.tbi',
            'merged_bincov': '.RD.txt.gz',
            'merged_bincov_index': '.RD.txt.gz.tbi',
            'SR_stats': '.SR.QC_matrix.txt',
            'PE_stats': '.PE.QC_matrix.txt',
            'BAF_stats': '.BAF.QC_matrix.txt',
            'RD_stats': '.RD.QC_matrix.txt',
            'median_cov': '_medianCov.transposed.bed',
            'merged_dels': '.DEL.bed.gz',
            'merged_dups': '.DUP.bed.gz',
            'Matrix_QC_plot': '.00_matrix_FC_QC.png',
        }
        for caller in SV_CALLERS:
            ending_by_key[f'std_{caller}_vcf_tar'] = f'.{caller}.tar.gz'

        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """Add jobs to Batch"""
        sids = dataset.get_sample_ids()

        input_by_sid = inputs.as_dict_by_target(stage=GatherSampleEvidence)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'ped_file': str(make_combined_ped(dataset)),
            'counts': [str(input_by_sid[sid]['coverage_counts']) for sid in sids],
            'SR_files': [str(input_by_sid[sid]['sr']) for sid in sids],
            'PE_files': [str(input_by_sid[sid]['pe']) for sid in sids],
            'SD_files': [str(input_by_sid[sid]['sd']) for sid in sids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(input_by_sid[sid][f'{caller}_vcf']) for sid in sids
            ]

        input_dict |= {
            'ref_copy_number_autosomal_contigs': 2,
            'allosomal_contigs': ['chrX', 'chrY'],
            'gcnv_qs_cutoff': 30,
            'min_svsize': 50,
            'run_matrix_qc': True,
            'matrix_qc_distance': 1000000,
        }
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
            ]
        )
        input_dict |= {
            'ref_dict': str(reference_path(f'broad/ref_fasta').with_suffix('.dict')),
        }

        # reference panel gCNV models
        input_dict |= get_ref_panel()

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
                'gcnv_gatk_docker',
                'cnmops_docker',
            ]
        )

        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


@stage(required_stages=GatherBatchEvidence)
class ClusterBatch(DatasetStage):
    """
    https://github.com/broadinstitute/gatk-sv#clusterbatch
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        * Clustered SV VCFs
        * Clustered depth-only call VCF
        * Metrics
        """

        ending_by_key = {
            'metrics_file_clusterbatch': '.metrics.tsv',
        }
        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'clustered_{caller}_vcf'] = f'.clustered-{caller}.vcf.gz'
            ending_by_key[
                f'clustered_{caller}_vcf_index'
            ] = f'.clustered-{caller}.vcf.gz.tbi'

        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Inputs:
        Standardized call VCFs (GatherBatchEvidence)
        Depth-only (DEL/DUP) calls (GatherBatchEvidence)
        """
        batch_evidence_d = inputs.as_dict(dataset, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'del_bed': str(batch_evidence_d[f'merged_dels']),
            'dup_bed': str(batch_evidence_d[f'merged_dups']),
            'ped_file': str(make_combined_ped(dataset)),
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcf_tar'] = str(
                batch_evidence_d[f'std_{caller}_vcf_tar']
            )

        input_dict |= {
            'depth_exclude_overlap_fraction': 0.5,
            'depth_interval_overlap': 0.8,
            'depth_clustering_algorithm': 'SINGLE_LINKAGE',
            'pesr_interval_overlap': 0.1,
            'pesr_breakend_window': 300,
            'pesr_clustering_algorithm': 'SINGLE_LINKAGE',
        }

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_pipeline_docker',
                'gatk_docker',
                'linux_docker',
                'sv_pipeline_base_docker',
            ]
        )

        input_dict |= get_references(
            [
                {'contig_list': 'primary_contigs_list'},
                {'depth_exclude_intervals': 'depth_exclude_list'},
                {'pesr_exclude_intervals': 'pesr_exclude_list'},
            ]
        )
        ref_fasta = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
        input_dict |= {
            'reference_fasta': str(ref_fasta),
            'reference_fasta_fai': str(ref_fasta) + '.fai',
            'reference_dict': str(ref_fasta.with_suffix('.dict')),
        }

        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


@stage(required_stages=[ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(DatasetStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Metrics file
        """

        ending_by_key = {
            'metrics': '.metrics.tsv',
            'metrics_common': '.metrics_common.tsv',
        }

        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        clusterbatch_d = inputs.as_dict(dataset, ClusterBatch)
        gatherbatchevidence_d = inputs.as_dict(dataset, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
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
            'ped_file': make_combined_ped(dataset),
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
            ]
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
                'rmsk',
                'segdups',
                {'autosome_contigs': 'autosome_file'},
                {'allosome_contigs': 'allosome_file'},
            ]
        )
        ref_fasta = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
        input_dict |= {
            'ref_dict': str(ref_fasta.with_suffix('.dict')),
        }

        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


@stage(required_stages=[GenerateBatchMetrics, ClusterBatch])
class FilterBatch(DatasetStage):
    """
    Filters poor quality variants and filters outlier samples. This workflow can
    be run all at once with the WDL at wdl/FilterBatch.wdl, or it can be run in three
    steps to enable tuning of outlier filtration cutoffs. The three sub-workflows are:

    * FilterBatchSites: Per-batch variant filtration
    * PlotSVCountsPerSample: Visualize SV counts per sample per type to help choose an
      IQR cutoff for outlier filtering, and preview outlier samples for a given cutoff
    * FilterBatchSamples: Per-batch outlier sample filtration; provide an appropriate
      outlier_cutoff_nIQR based on the SV count plots and outlier previews from step 2.
      Note that not removing high outliers can result in increased compute cost and
      a higher false positive rate in later steps.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        * Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        * Filtered depth-only call VCF with outlier samples excluded
        * Random forest cutoffs file
        * PED file with outlier samples excluded"""
        ending_by_key = {
            'done': '.done',
        }

        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        metrics_d = inputs.as_dict(dataset, GenerateBatchMetrics)
        clusterbatch_d = inputs.as_dict(dataset, ClusterBatch)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'ped_file': make_combined_ped(dataset),
            'evidence_metrics': metrics_d['metrics'],
            'evidence_metrics_common': metrics_d['metrics_common'],
            'outlier_cutoff_nIQR': '6',
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_base_docker',
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'linux_docker',
            ]
        )

        input_dict |= get_references(
            [
                'primary_contigs_list',
            ]
        )

        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


@stage(required_stages=[FilterBatch])
class GenotypeBatch(DatasetStage):
    """
    Genotypes a batch of samples across unfiltered variants combined across all batches.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """
        return {}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        pass


@stage(required_stages=GenotypeBatch)
class MakeCohortVcf(DatasetStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        pass


@stage(required_stages=MakeCohortVcf)
class AnnotateVcf(DatasetStage):
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

    def expected_outputs(self, dataset: Dataset) -> dict:
        vcf_path = str(self.tmp_prefix / 'sv' / f'{dataset.name}.vcf.gz')
        return {
            'vcf': vcf_path,
            'tbi': vcf_path + '.tbi',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        pass
