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
from cpg_utils.config import get_config, ConfigError
from cpg_utils.hail_batch import command, reference_path, image_path
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


GATK_SV_COMMIT = 'a73237cf9d9e321df3aa81c890def7b504a25c7f'
SV_CALLERS = ['manta', 'wham', 'scramble']
_FASTA = None


def get_fasta() -> Path:
    """
    find or return the fasta to use

    Returns:

    """
    global _FASTA
    if _FASTA is None:
        _FASTA = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
    return _FASTA


def get_images(keys: list[str], allow_missing=False) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.

    Args:
        keys (list): all the images to get
        allow_missing (bool): if False, require all query keys to be found

    Returns:
        dict of image keys to image paths
        or AssertionError
    """

    if not allow_missing:
        image_keys = get_config()['images'].keys()
        query_keys = set(keys)
        if not query_keys.issubset(image_keys):
            raise KeyError(f'Unknown image keys: {query_keys - image_keys}')

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
        except ConfigError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))

    return res


def add_gatk_sv_jobs(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path | list[Path]],
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
    for key, value in expected_out_dict.items():
        if isinstance(value, list):
            outputs_to_collect[key] = CromwellOutputType.array_path(
                name=f'{wfl_name}.{key}', length=len(value)
            )
        else:
            outputs_to_collect[key] = CromwellOutputType.single_path(
                f'{wfl_name}.{key}'
            )

    driver_image = driver_image or image_path('cpg_workflows')

    # pre-process input_dict
    paths_as_strings: dict = {}
    for key, value in input_dict.items():
        if isinstance(value, Path):
            paths_as_strings[f'{wfl_name}.{key}'] = str(value)
        elif isinstance(value, (list, set)):
            paths_as_strings[f'{wfl_name}.{key}'] = [str(v) for v in value]
        else:
            paths_as_strings[f'{wfl_name}.{key}'] = value

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
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(driver_image)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        if isinstance(resource, list):
            for source, dest in zip(resource, out_path):
                cmds.append(f'gsutil cp "$(cat {source})" "{dest}"')
        else:
            cmds.append(f'gsutil cp "$(cat {resource})" "{out_path}"')
    copy_j.command(command(cmds, setup_gcp=True))
    return [submit_j, copy_j]


def get_ref_panel(keys: list[str] | None = None) -> dict:
    return {
        k: v
        for k, v in {
            'ref_panel_samples': get_config()['sv_ref_panel']['ref_panel_samples'],
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
        ending_by_key['pesr_split'] = 'sr.txt.gz'
        ending_by_key['pesr_split_index'] = 'sr.txt.gz.tbi'
        # discordant paired end reads:
        ending_by_key['pesr_disc'] = 'pe.txt.gz'
        ending_by_key['pesr_disc_index'] = 'pe.txt.gz.tbi'
        # site depth:
        ending_by_key['pesr_sd'] = 'sd.txt.gz'
        ending_by_key['pesr_sd_index'] = 'sd.txt.gz.tbi'

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
            'reference_fasta': str(get_fasta()),
            'reference_index': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'reference_version': '38'
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
            'bincov_median': f'{dataset.name}_medianCov.transposed.bed'
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
            'SR_files': [str(input_by_sid[sid]['pesr_split']) for sid in sids],
            'PE_files': [str(input_by_sid[sid]['pesr_disc']) for sid in sids],
            'SD_files': [str(input_by_sid[sid]['pesr_sd']) for sid in sids],
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
                str(input_by_sid[sid][f'{caller}_vcf']) for sid in sids
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
            ]
        )

        # reference panel gCNV models
        input_dict |= get_ref_panel()

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'linux_docker',
                'condense_counts_docker',
                'gatk_docker',
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
            'del_bed': str(batch_evidence_d['merged_dels']),
            'dup_bed': str(batch_evidence_d['merged_dups']),
            'ped_file': str(make_combined_ped(dataset)),
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
            input_dict[f'{caller}_vcf_tar'] = str(
                batch_evidence_d[f'std_{caller}_vcf_tar']
            )

        input_dict |= get_images(
            ['sv_base_mini_docker', 'sv_pipeline_docker', 'gatk_docker', 'linux_docker']
        )

        input_dict |= get_references(
            [
                {'contig_list': 'primary_contigs_list'},
                {'depth_exclude_intervals': 'depth_exclude_list'},
                {'pesr_exclude_intervals': 'pesr_exclude_list'},
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


@stage(required_stages=[ClusterBatch, GatherBatchEvidence])
class GenerateBatchMetrics(DatasetStage):
    """
    Generates variant metrics for filtering.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Metrics files
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
            'ref_dict': str(get_fasta().with_suffix('.dict')),
        }

        for caller in SV_CALLERS + ['depth']:
            input_dict[f'{caller}_vcf'] = clusterbatch_d[f'clustered_{caller}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'sv_base_docker',
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
        ending_by_key: dict = {
            'metrics_file_filterbatch': '.metrics.tsv',
            'filtered_pesr_vcf': '.filtered-pesr-merged.vcf.gz',
            'cutoffs': '.cutoffs',
            'scores': '.updated_scores',
            'RF_intermediate_files': '.RF_intermediate_files.tar.gz',
            'outlier_samples_excluded_file': '.outliers.samples.list',
            'batch_samples_postOutlierExclusion_file': '.outliers_excluded.samples.list',
        }
        for caller in SV_CALLERS + ['depth']:
            ending_by_key[f'filtered_{caller}_vcf'] = f'.filtered-{caller}.vcf.gz'

            # unsure why, scramble doesn't export this file
            if caller != 'scramble':
                ending_by_key[
                    f'sites_filtered_{caller}_vcf'
                ] = f'.sites-filtered-{caller}.vcf.gz'

        ending_by_key['sv_counts'] = [
            f'.{caller}.with_evidence.svcounts.txt' for caller in SV_CALLERS + ['depth']
        ]
        ending_by_key['sv_count_plots'] = [
            f'.{caller}.with_evidence.all_SVTYPEs.counts_per_sample.png'
            for caller in SV_CALLERS + ['depth']
        ]
        d: dict[str, Path | list[Path]] = {}
        for key, ending in ending_by_key.items():
            if isinstance(ending, str):
                fname = f'{dataset.name}{ending}'
                d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
            elif isinstance(ending, list):
                d[key] = [
                    dataset.prefix()
                    / 'gatk_sv'
                    / self.name.lower()
                    / f'{dataset.name}{e}'
                    for e in ending
                ]
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


@stage(required_stages=[FilterBatch, GatherBatchEvidence])
class GenotypeBatch(DatasetStage):
    """
    Genotypes a batch of samples across filtered variants combined across all batches.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
        Filtered depth-only call VCF with outlier samples excluded
        PED file with outlier samples excluded
        List of SR pass variants
        List of SR fail variants
        """
        ending_by_key = {
            'sr_bothside_pass': '.genotype_SR_part2_bothside_pass.txt',
            'sr_background_fail': '.genotype_SR_part2_background_fail.txt',
            'trained_PE_metrics': '.pe_metric_file.txt',
            'trained_SR_metrics': '.sr_metric_file.txt',
            'regeno_coverage_medians': '.regeno.coverage_medians_merged.bed',
            'metrics_file_genotypebatch': '.metrics.tsv',
        }

        # really unsure about this bit - it looks like both inputs are
        # the same prior output? TrainRDGenotyping.pesr/depth_sepcutoff
        for mode in ['pesr', 'depth']:
            ending_by_key |= {
                f'trained_genotype_{mode}_pesr_sepcutoff': f'.{mode}.pesr_sepcutoff.txt',
                f'trained_genotype_{mode}_depth_sepcutoff': f'.{mode}.depth_sepcutoff.txt',
                f'genotyped_{mode}_vcf': f'.{mode}.vcf.gz',
                f'genotyped_{mode}_vcf_index': f'.{mode}.vcf.gz.tbi',
            }
        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:

        filterbatch_d = inputs.as_dict(dataset, FilterBatch)
        batchevidence_d = inputs.as_dict(dataset, GatherBatchEvidence)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'ped_file': make_combined_ped(dataset),
            'n_per_split': 5000,
            'n_RD_genotype_bins': 100000,
            'coveragefile': batchevidence_d['merged_bincov'],  # unsure
            'coveragefile_index': batchevidence_d['merged_bincov_index'],  # unsure
            'discfile': batchevidence_d['merged_PE'],
            'discfile_index': batchevidence_d['merged_PE_index'],
            'splitfile': batchevidence_d['merged_SR'],
            'splitfile_index': batchevidence_d['merged_SR_index'],
            'medianfile': batchevidence_d['median_cov'],
            'rf_cutoffs': filterbatch_d['cutoffs'],
            # instead of pulling from reference:
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'reference_build': 'hg38'
        }

        for mode in ['pesr', 'depth']:
            input_dict[f'batch_{mode}_vcf'] = filterbatch_d[f'filtered_{mode}_vcf']
            input_dict[f'cohort_{mode}_vcf'] = filterbatch_d[f'filtered_{mode}_vcf']

        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'linux_docker',
            ]
        )
        input_dict |= get_references(
            [
                'primary_contigs_list',
                'bin_exclude',
                'seed_cutoffs',
                'pesr_exclude_list'
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


# @stage(required_stages=GenotypeBatch)
# class RegenotypeCNVs(DatasetStage):
#     """
#     Optional, Re-genotypes probable mosaic variants across multiple batches.
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict:
#         pass
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
#         pass


@stage(required_stages=[FilterBatch, GenotypeBatch, GatherBatchEvidence])
class MakeCohortVcf(DatasetStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    If RegenotypeCNVs is run, this stage will use the output of that stage as input.
    Otherwise, it will use the output of GenotypeBatch.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """

        Args:
            dataset ():

        Returns:

        """
        ending_by_key = {
            'vcf': '.cleaned.vcf.gz',
            'vcf_index': '.cleaned.vcf.gz.tbi',
            'vcf_qc': '.cleaned_SV_VCF_QC_output.tar.gz',
            # if merge_intermediate_vcfs is enabled
            # 'cluster_vcf': '.combine_batches.vcf.gz',
            # 'cluster_vcf_index': '.combine_batches.vcf.gz.tbi',
            # 'complex_resolve_vcf': '.complex_resolve.vcf.gz',
            # 'complex_resolve_vcf_index': '.complex_resolve.vcf.gz.tbi',
            # 'complex_genotype_vcf': '.complex_genotype.vcf.gz',
            # 'complex_genotype_vcf_index': '.complex_genotype.vcf.gz.tbi',
            'metrics_file_makecohortvcf': '.metrics.tsv'
        }
        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """

        Args:
            dataset ():
            inputs ():

        Returns:

        """
        batchevidence_d = inputs.as_dict(dataset, GatherBatchEvidence)
        filterbatch_d = inputs.as_dict(dataset, FilterBatch)
        genotypebatch_d = inputs.as_dict(dataset, GenotypeBatch)

        input_dict: dict[str, Any] = {
            'cohort_name': dataset.name,
            'batches': [dataset.name],
            'ped_file': make_combined_ped(dataset),
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
            'pesr_vcfs': [genotypebatch_d['genotyped_pesr_vcf']],
            'depth_vcfs': [genotypebatch_d['genotyped_depth_vcf']],
            'disc_files': [batchevidence_d['merged_PE']],
            'bincov_files': [batchevidence_d['merged_bincov']],
            'raw_sr_bothside_pass_files': [genotypebatch_d['sr_bothside_pass']],
            'raw_sr_background_fail_files': [genotypebatch_d['sr_background_fail']],
            'depth_gt_rd_sep_files': [genotypebatch_d['trained_genotype_depth_depth_sepcutoff']],
            'median_coverage_files': [batchevidence_d['median_cov']],
            'rf_cutoff_files': [filterbatch_d['cutoffs']],
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
                'sv_base_mini_docker',
                'linux_docker',
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
        # vcf_path = str(self.tmp_prefix / 'sv' / f'{dataset.name}.vcf.gz')
        ending_by_key = {
            'output_vcf': '.annotated.vcf.gz',
            'output_vcf_idx': '.annotated.vcf.gz.tbi',
        }
        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:

        make_vcf_d = inputs.as_dict(dataset, MakeCohortVcf)

        input_dict: dict[str, Any] = {
            'prefix': dataset.name,
            'vcf': make_vcf_d['vcf'],
            'vcf_idx': make_vcf_d['vcf_index'],
            'ped_file': make_combined_ped(dataset),
            'sv_per_shard': 5000,
            'max_shards_per_chrom_step1': 200,
            'min_records_per_shard_step1': 5000
        }
        input_dict |= get_references(
            [
                'protein_coding_gtf',
                {'contig_list': 'primary_contigs_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'gatk_docker',
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
