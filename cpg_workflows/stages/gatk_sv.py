"""
Driver for computing structural variants from GATK-SV from WGS data.
"""
from os.path import join
from typing import Any

from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch import Resource
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
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

GATK_SV_COMMIT = '50d47e68ad05b29c78cba16606974ef0dca1b113'
SV_CALLERS = ['manta', 'wham', 'scramble']
REF_FASTA_KEY = 'broad/bwa_ref_fasta'  # change to 'broad/ref_fasta' for DRAGMAP CRAMs


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
            res[key] = str(reference_path(f'broad/sv/resources/{ref_d_key}'))
        except KeyError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))

    return res


def add_gatk_sv_job(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path],
    sample_id: str | None = None,
) -> Job:
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

    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset.name)
    j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
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
        input_dict={f'{wfl_name}.{k}': v for k, v in input_dict.items()},
        outputs_to_collect=outputs_to_collect,
        driver_image=image_path('hail'),
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(image_path('hail'))
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        cmds.append(f'gsutil cp $(cat {resource}) {out_path}')
    copy_j.command(command(cmds, setup_gcp=True))
    return copy_j


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
        input_dict |= {
            'reference_fasta': str(ref_fasta := reference_path(REF_FASTA_KEY)),
            'reference_index': str(ref_fasta) + '.fai',
            'reference_dict': str(ref_fasta.with_suffix('.dict')),
        }
        input_dict['reference_version'] = '38'

        expected_d = self.expected_outputs(sample)

        j = add_gatk_sv_job(
            batch=get_batch(),
            dataset=sample.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
        )
        return self.make_outputs(sample, data=expected_d, jobs=[j])


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
        """Add jobs to batch"""
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
        j = add_gatk_sv_job(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=[j])


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

        # PED file must include reference panel samples as well, so concatenating
        # the `dataset` PED file with the reference panel PED file:
        combined_ped_path = (
            self.tmp_prefix
            / 'ped_with_ref_panel'
            / f'{dataset.alignment_inputs_hash()}.ped'
        )
        with combined_ped_path.open('w') as out:
            with dataset.write_ped_file().open() as f:
                out.write(f.read())
            # The ref panel PED doesn't have a header, so can safely concatenate:
            with reference_path('broad/sv/ref_panel/ped_file').open() as f:
                out.write(f.read())

        input_by_sid = inputs.as_dict_by_target(stage=GatherSampleEvidence)

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'ped_file': str(combined_ped_path),
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
        input_dict |= {
            'ref_panel_samples': get_config()['sv_ref_panel']['ref_panel_samples'],
            'ref_panel_bincov_matrix': str(
                reference_path('broad/sv/ref_panel/ref_panel_bincov_matrix')
            ),
            'contig_ploidy_model_tar': str(
                reference_path('broad/sv/ref_panel/contig_ploidy_model_tar')
            ),
            'gcnv_model_tars': [
                str(reference_path('broad/sv/ref_panel/model_tar_tmpl')).format(shard=i)
                for i in range(get_config()['sv_ref_panel']['model_tar_cnt'])
            ],
            'ref_panel_PE_files': [
                str(reference_path('broad/sv/ref_panel/ref_panel_PE_file_tmpl')).format(
                    sample=s
                )
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SR_files': [
                str(reference_path('broad/sv/ref_panel/ref_panel_SR_file_tmpl')).format(
                    sample=s
                )
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SD_files': [
                str(reference_path('broad/sv/ref_panel/ref_panel_SD_file_tmpl')).format(
                    sample=s
                )
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
        }

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
        j = add_gatk_sv_job(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=[j])
