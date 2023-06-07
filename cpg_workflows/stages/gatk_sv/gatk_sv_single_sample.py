"""
single-sample components of the GATK SV workflow
"""
import logging
from typing import Any

from cpg_utils import Path
from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import (
    stage,
    SampleStage,
    StageOutput,
    StageInput,
    Sample,
    Cohort,
    CohortStage,
)
from cpg_workflows.jobs import sample_batching

from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
    _sv_individual_meta,
    SV_CALLERS,
)
from cpg_utils.config import get_config


@stage(analysis_type='sv', update_analysis_meta=_sv_individual_meta)
class GatherSampleEvidence(SampleStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Expected to produce coverage counts, a VCF for each variant caller,
        and a txt for each type of SV evidence (SR, PE, SD).

        it's critical to separate the ending with a dot, e.g.: `*.sr.txt.gz`,
        These files are passed to `gatk PrintSVEvidence`, that determines file
        format based on the file name.
        It would strongly expect the files to end exactly with either
        `.sr.txt.gz`, `.pe.txt.gz`, or `.sd.txt.gz`, otherwise it would fail with
        "A USER ERROR has occurred: Cannot read file:///cromwell_root/... because
        no suitable codecs found".
        """
        d: dict[str, Path] = {
            'coverage_counts': sample.make_sv_evidence_path
            / f'{sample.id}.coverage_counts.tsv.gz',
            # split reads
            'pesr_split': sample.make_sv_evidence_path / f'{sample.id}.sr.txt.gz',
            'pesr_split_index': sample.make_sv_evidence_path
            / f'{sample.id}.sr.txt.gz.tbi',
            # site depth
            'pesr_sd': sample.make_sv_evidence_path / f'{sample.id}.sd.txt.gz',
            'pesr_sd_index': sample.make_sv_evidence_path
            / f'{sample.id}.sd.txt.gz.tbi',
            # discordant paired reads
            'pesr_disc': sample.make_sv_evidence_path / f'{sample.id}.pe.txt.gz',
            'pesr_disc_index': sample.make_sv_evidence_path
            / f'{sample.id}.pe.txt.gz.tbi',
        }

        # Caller's VCFs
        for caller in SV_CALLERS:
            d[f'{caller}_vcf'] = (
                sample.make_sv_evidence_path / f'{sample.id}.{caller}.vcf.gz'
            )
            d[f'{caller}_index'] = (
                sample.make_sv_evidence_path / f'{sample.id}.{caller}.vcf.gz.tbi'
            )

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
            'reference_version': '38',
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
class EvidenceQC(CohortStage):
    """
    https://github.com/broadinstitute/gatk-sv#evidenceqc
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        fname_by_key = {
            'ploidy_matrix': 'ploidy_matrix.bed.gz',
            'ploidy_plots': 'ploidy_plots.tar.gz',
            'WGD_dist': 'WGD_score_distributions.pdf',
            'WGD_matrix': 'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores': 'WGD_scores.txt.gz',
            'bincov_matrix': 'RD.txt.gz',
            'bincov_matrix_index': 'RD.txt.gz.tbi',
            'bincov_median': 'medianCov.transposed.bed',
            'qc_table': 'evidence_qc_table.tsv',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        return {key: self.prefix / fname for key, fname in fname_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        d = inputs.as_dict_by_target(GatherSampleEvidence)
        sids = cohort.get_sample_ids()

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'samples': sids,
            'run_vcf_qc': True,
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

        input_dict |= get_references(['genome_file', 'wgd_scoring_mask'])

        expected_d = self.expected_outputs(cohort)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(cohort, data=expected_d, jobs=jobs)


@stage(required_stages=EvidenceQC, analysis_type='sv', analysis_keys=['batch_json'])
class CreateSampleBatches(CohortStage):
    """
    uses the values generated in EvidenceQC, does some clustering
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {'batch_json': self.prefix / '{pcr_status}_batches.json'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        expected = self.expected_outputs(cohort)
        pcr_plus, pcr_neg = [], []
        for sample in cohort.get_samples():
            if sample.meta.get('pcr_status', 'unknown') == 'PCR-':
                pcr_neg.append(sample.id)
            else:
                pcr_plus.append(sample.id)

        # get those settings
        min_batch_size = get_config()['workflow'].get('min_batch_size', 100)
        max_batch_size = get_config()['workflow'].get('max_batch_size', 300)

        for status, samples in [('PCRMINUS', pcr_neg), ('PCRPLUS', pcr_plus)]:
            if len(samples) < min_batch_size:
                logging.info(f'Too few {status} samples to form batches')
                return

            # I think this can just be a PythonJob?
            py_job = get_batch().new_python_job(f'create_{status}_sample_batches')
            py_job.image(get_config()['workflow']['driver_image'])
            py_job.call(
                sample_batching.partition_batches,
                inputs.as_dict(cohort, EvidenceQC)['qc_table'],
                samples,
                expected['batch_json'],
                min_batch_size,
                max_batch_size,
            )

        return self.make_outputs(cohort, data=expected, jobs=py_job)
