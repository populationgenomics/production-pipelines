"""
Driver for computing structural variants from GATK-SV from WGS data.
"""
from typing import Any

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, genome_build
from cpg_workflows.batch import get_batch
from cpg_workflows.stages.gatk_sv import (
    GatherSampleEvidence,
    get_ref_panel,
    get_images,
    get_references,
    add_gatk_sv_jobs,
)
from cpg_workflows.workflow import (
    stage,
    SampleStage,
    StageOutput,
    StageInput,
    Sample,
)

GATK_SV_COMMIT = '7444972aabf9fbb19ac3e5cfc9387098f5a411d2'
SV_CALLERS = ['manta', 'wham', 'scramble']


@stage(required_stages=GatherSampleEvidence)
class GATKSVPipelineSingleSample(SampleStage):
    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        d: dict[str, Path] = {}

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
        assert sample.cram, sample

        d = inputs.as_dict_by_target(GatherSampleEvidence)

        input_dict: dict[str, Any] = {
            'case_counts_file': str(d['coverage_counts']),
        }
        for caller in SV_CALLERS:
            input_dict[f'case_{caller}_vcf'] = str(d[f'{caller}_vcf'])

        input_dict |= {
            'batch': sample.dataset.name,
            'sample': sample.id,
            # GatherSampleEvidence
            'reference_version': genome_build()[:-2],
            # EvidenceQC
            'run_vcf_qc': True,  # generates <caller>_qc_low/<caller>_qc_high
            # GatherBatchEvidence
            'min_svsize': 50,
            'ref_copy_number_autosomal_contigs': 2,
            'allosomal_contigs': ['chrX', 'chrY'],
            'gcnv_qs_cutoff': 30,
            'matrix_qc_distance': 1000000,
            # ClusterBatch
        }
        # reference panel gCNV models
        input_dict |= get_ref_panel(
            [
                'ref_ped_file',
                'ref_samples_list',
                'ref_panel_vcf',
                'contig_ploidy_model_tar',
                'gcnv_model_tars',
                'ref_panel_bincov_matrix',
                'cnmops_exclude_list',
            ]
        )
        with (_path := self.tmp_prefix / 'gcnv_model_tars_list.txt').open('w') as f:
            for _tar_path in input_dict['gcnv_model_tars']:
                f.write(f'{_tar_path}\n')
        input_dict['gcnv_model_tars_list'] = str(_path)
        del input_dict['gcnv_model_tars']

        with (_path := self.tmp_prefix / 'ref_pesr_disc_files_list.txt').open('w') as f:
            for _tar_path in input_dict['ref_panel_PE_files']:
                f.write(f'{_tar_path}\n')
        input_dict['ref_pesr_disc_files_list'] = str(_path)
        del input_dict['ref_panel_PE_files']

        with (_path := self.tmp_prefix / 'ref_pesr_split_files_list.txt').open(
            'w'
        ) as f:
            for _tar_path in input_dict['ref_panel_SR_files']:
                f.write(f'{_tar_path}\n')
        input_dict['ref_pesr_split_files_list'] = str(_path)
        del input_dict['ref_panel_SR_files']

        with (_path := self.tmp_prefix / 'ref_pesr_sd_files_list.txt').open('w') as f:
            for _tar_path in input_dict['ref_panel_SD_files']:
                f.write(f'{_tar_path}\n')
        input_dict['ref_pesr_sd_files_list'] = str(_path)
        del input_dict['ref_panel_SD_files']

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_hail_docker',
                'sv_pipeline_updates_docker',
                'sv_pipeline_rdtest_docker',
                'sv_pipeline_base_docker',
                'sv_pipeline_qc_docker',
                'linux_docker',
                'cnmops_docker',
                'gatk_docker',
                'condense_counts_docker',
                'genomes_in_the_cloud_docker',
                'samtools_cloud_docker',
                'cloud_sdk_docker',
            ]
        )
        input_dict |= get_references(
            [
                'genome_file',
                'primary_contigs_list',
                'primary_contigs_fai',
                'autosome_file',
                'allosome_file',
                # GatherSampleEvidence (not used because already ran)
                'preprocessed_intervals',
                'manta_region_bed',
                {'sd_locs_vcf': 'dbsnp_vcf'},
                'wham_include_list_bed_file',
                # EvidenceQC
                'wgd_scoring_mask',
                # GatherBatchEvidence
                'contig_ploidy_model_tar',
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
