"""
single-sample components of the GATK SV workflow
"""

import logging
from typing import Any

from cpg_utils import Path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, try_get_ar_guid
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import sample_batching
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    SV_CALLERS,
    CromwellJobSizes,
    _sv_individual_meta,
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
)
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


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
            'coverage_counts': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.coverage_counts.tsv.gz',
            # split reads
            'pesr_split': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sr.txt.gz',
            'pesr_split_index': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sr.txt.gz.tbi',
            # site depth
            'pesr_sd': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sd.txt.gz',
            'pesr_sd_index': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.sd.txt.gz.tbi',
            # discordant paired reads
            'pesr_disc': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.pe.txt.gz',
            'pesr_disc_index': sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.pe.txt.gz.tbi',
        }

        # Caller's VCFs
        for caller in SV_CALLERS:
            d[f'{caller}_vcf'] = sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.{caller}.vcf.gz'
            d[f'{caller}_index'] = sequencing_group.make_sv_evidence_path / f'{sequencing_group.id}.{caller}.vcf.gz.tbi'

        return d

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        Add jobs to batch
        Adds billing-related labels to the Cromwell job(s)
        """
        assert sequencing_group.cram, sequencing_group

        input_dict: dict[str, Any] = dict(
            bam_or_cram_file=str(sequencing_group.cram),
            bam_or_cram_index=str(sequencing_group.cram) + '.crai',
            sample_id=sequencing_group.id,
            reference_fasta=str(get_fasta()),
            reference_index=str(get_fasta()) + '.fai',
            reference_dict=str(get_fasta().with_suffix('.dict')),
            reference_version='38',
        )

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
            ],
        )
        input_dict |= get_references(
            [
                'primary_contigs_fai',
                'primary_contigs_list',
                'preprocessed_intervals',
                'manta_region_bed',
                'wham_include_list_bed_file',
                {'sd_locs_vcf': 'dbsnp_vcf'},
            ],
        )

        expected_d = self.expected_outputs(sequencing_group)

        # billing labels!
        # https://cromwell.readthedocs.io/en/stable/wf_options/Google/
        # these must conform to the regex [a-z]([-a-z0-9]*[a-z0-9])?
        billing_labels = {
            'dataset': sequencing_group.dataset.name,  # already lowercase
            'sequencing-group': sequencing_group.id.lower(),
            'stage': self.name.lower(),
            AR_GUID_NAME: try_get_ar_guid(),
        }

        # add some max-polling interval jitter for each sample
        # cromwell_status_poll_interval is a number (int, seconds)
        # this is used to determine how often to poll Cromwell for completion status
        # we alter the per-sample maximum to be between 5 and 30 minutes for this
        # long-running job, so samples poll on different intervals, spreading load
        jobs = add_gatk_sv_jobs(
            dataset=sequencing_group.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sequencing_group_id=sequencing_group.id,
            labels=billing_labels,
            job_size=CromwellJobSizes.LARGE,
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
            input_dict[f'{caller}_vcfs'] = [str(d[sid][f'{caller}_vcf']) for sid in sgids]

        input_dict |= get_images(
            ['sv_base_mini_docker', 'sv_base_docker', 'sv_pipeline_docker', 'sv_pipeline_qc_docker'],
        )

        input_dict |= get_references(['genome_file', 'wgd_scoring_mask'])

        expected_d = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        # runs for approx 5 hours, depending on sample count
        jobs = add_gatk_sv_jobs(
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            labels=billing_labels,
            job_size=CromwellJobSizes.MEDIUM,
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
        return {'batch_json': self.analysis_prefix / 'sgid_batches.json'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        this stage has been clarified - this will only run on Genomes
        Exomes were never a supported case, they have a separate pipeline
        """

        expected = self.expected_outputs(cohort)

        if config_retrieve(['workflow', 'sequencing_type']) != 'genome':
            raise RuntimeError('This workflow is not intended for Exome data')

        sequencing_groups = [sequencing_group.id for sequencing_group in cohort.get_sequencing_groups()]

        # get those settings
        min_batch_size = config_retrieve(['workflow', 'min_batch_size'], 100)
        max_batch_size = config_retrieve(['workflow', 'max_batch_size'], 300)

        if len(sequencing_groups) < min_batch_size:
            logging.error('Too few sequencing groups to form batches')
            raise RuntimeError('too few samples to create batches')

        py_job = get_batch().new_python_job('create_sample_batches')
        py_job.image(config_retrieve(['workflow', 'driver_image']))
        py_job.call(
            sample_batching.partition_batches,
            inputs.as_dict(cohort, EvidenceQC)['qc_table'],
            sequencing_groups,
            expected['batch_json'],
            min_batch_size,
            max_batch_size,
        )

        return self.make_outputs(cohort, data=expected, jobs=py_job)
