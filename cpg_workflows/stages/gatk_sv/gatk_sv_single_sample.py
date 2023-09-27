"""
single-sample components of the GATK SV workflow
"""
import logging
from typing import Any

from cpg_utils import Path
from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import (
    stage,
    SequencingGroupStage,
    StageOutput,
    StageInput,
    SequencingGroup,
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


@stage(
    analysis_keys=[f'{caller}_vcf' for caller in SV_CALLERS],
    analysis_type='sv',
    update_analysis_meta=_sv_individual_meta,
)
class GatherSampleEvidence(SequencingGroupStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence
    https://github.com/broadinstitute/gatk-sv/blob/master/wdl/GatherSampleEvidence.wdl
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
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
            'coverage_counts': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.coverage_counts.tsv.gz',
            # split reads
            'pesr_split': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.sr.txt.gz',
            'pesr_split_index': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.sr.txt.gz.tbi',
            # site depth
            'pesr_sd': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.sd.txt.gz',
            'pesr_sd_index': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.sd.txt.gz.tbi',
            # discordant paired reads
            'pesr_disc': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.pe.txt.gz',
            'pesr_disc_index': sequencing_group.make_sv_evidence_path
            / f'{sequencing_group.id}.pe.txt.gz.tbi',
        }

        # Caller's VCFs
        for caller in SV_CALLERS:
            d[f'{caller}_vcf'] = (
                sequencing_group.make_sv_evidence_path
                / f'{sequencing_group.id}.{caller}.vcf.gz'
            )
            d[f'{caller}_index'] = (
                sequencing_group.make_sv_evidence_path
                / f'{sequencing_group.id}.{caller}.vcf.gz.tbi'
            )

        return d

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput:
        """
        Add jobs to batch
        Adds billing-related labels to the Cromwell job(s)
        """
        assert sequencing_group.cram, sequencing_group

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(sequencing_group.cram),
            'bam_or_cram_index': str(sequencing_group.cram) + '.crai',
            'sample_id': sequencing_group.id,
            'reference_fasta': str(get_fasta()),
            'reference_index': str(get_fasta()) + '.fai',
            'reference_dict': str(get_fasta().with_suffix('.dict')),
            'reference_version': '38',
        }

        input_dict |= get_images(
            [
                'sv_pipeline_base_docker',
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'samtools_cloud_docker',
                'manta_docker',
                'scramble_docker',
                'wham_docker',
                'gatk_docker',
                'gatk_docker_pesr_override',
                'genomes_in_the_cloud_docker',
                'cloud_sdk_docker',
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

        # find any additional arguments to pass to Cromwell
        if override := get_config()['resource_overrides'].get('GatherSampleEvidence'):
            input_dict |= override

        expected_d = self.expected_outputs(sequencing_group)

        # billing labels!
        billing_labels = {
            'dataset': sequencing_group.dataset.name,
            'sequencing_group': sequencing_group.id,
        }

        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=sequencing_group.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sequencing_group_id=sequencing_group.id,
            labels=billing_labels,
        )
        return self.make_outputs(sequencing_group, data=expected_d, jobs=jobs)


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
            'bincov_matrix': f'{self.name}.RD.txt.gz',
            'bincov_matrix_index': f'{self.name}.RD.txt.gz.tbi',
            'bincov_median': 'medianCov.transposed.bed',
            'qc_table': 'evidence_qc_table.tsv',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        return {key: self.prefix / fname for key, fname in fname_by_key.items()}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        d = inputs.as_dict_by_target(GatherSampleEvidence)
        sgids = cohort.get_sequencing_group_ids()

        input_dict: dict[str, Any] = {
            'batch': cohort.name,
            'samples': sgids,
            'run_vcf_qc': True,
            'counts': [str(d[sid]['coverage_counts']) for sid in sgids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(d[sid][f'{caller}_vcf']) for sid in sgids
            ]

        input_dict |= get_images(
            [
                'sv_base_mini_docker',
                'sv_base_docker',
                'sv_pipeline_docker',
                'sv_pipeline_qc_docker',
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
    uses the values generated in EvidenceQC
    splits the sequencing groups into batches based on median coverage,
    PCR +/- status, and Sex

    The output of this Stage will contain the distinct SG batches to use for the
    following series of Stages. For now, standard practice is to create a separate
    minimal configuration file for each sub-batch, containing the list of SG IDs
    as the `only_sgs` key. The gatk_sv_multisample_1 and gatk_sv_sandwich WFs are
    then run separately for each sub-batch, with the active SGs controlled via the
    config contents.

    When we move to custom cohorts, the output of this stage will be used as input
    when generating a custom Metamist cohort per sub-batch.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {'batch_json': self.prefix / 'pcr_{pcr_status}_batches.json'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        expected = self.expected_outputs(cohort)

        # PCR +/- logic is only relevant to exomes
        if get_config()['workflow'].get('sequencing_type') != 'genome':
            sequencing_group_types = {
                'exome': [sg.id for sg in cohort.get_sequencing_groups()]
            }
        # within exomes, divide into PCR- and all other sequencing groups
        else:
            sequencing_group_types = {'positive': [], 'negative': []}
            for sequencing_group in cohort.get_sequencing_groups():
                if sequencing_group.meta.get('pcr_status', 'unknown') == 'negative':
                    sequencing_group_types['negative'].append(sequencing_group.id)
                else:
                    sequencing_group_types['positive'].append(sequencing_group.id)

        # get those settings
        min_batch_size = get_config()['workflow'].get('min_batch_size', 100)
        max_batch_size = get_config()['workflow'].get('max_batch_size', 300)

        all_jobs = []

        for status, sequencing_groups in sequencing_group_types.items():
            if not sequencing_groups:
                logging.info(f'No {status} sequencing groups found')
                continue
            elif len(sequencing_groups) < min_batch_size:
                logging.info(f'Too few {status} sequencing groups to form batches')
                continue
            logging.info(f'Creating {status} sequencing group batches')
            py_job = get_batch().new_python_job(f'create_{status}_sample_batches')
            py_job.image(get_config()['workflow']['driver_image'])
            py_job.call(
                sample_batching.partition_batches,
                inputs.as_dict(cohort, EvidenceQC)['qc_table'],
                sequencing_groups,
                expected['batch_json'],
                min_batch_size,
                max_batch_size,
            )
            all_jobs.append(py_job)

        return self.make_outputs(cohort, data=expected, jobs=all_jobs)
