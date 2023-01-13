"""
Metamist wrapper to get input samples.
"""

import logging

from cpg_utils.config import get_config, update_dict

from .metamist import get_metamist, Sequence, AnalysisType, MetamistError
from .targets import Cohort, Sex, PedigreeInfo


_cohort: Cohort | None = None


def get_cohort() -> Cohort:
    """Return the cohort object"""
    global _cohort
    if not _cohort:
        _cohort = create_cohort()
    return _cohort


def create_cohort() -> Cohort:
    """
    Add datasets in the cohort. There exists only one cohort for the workflow run.
    """
    analysis_dataset_name = get_config()['workflow']['dataset']
    dataset_names = get_config()['workflow'].get(
        'input_datasets', [analysis_dataset_name]
    )
    skip_datasets = get_config()['workflow'].get('skip_datasets', [])
    dataset_names = [d for d in dataset_names if d not in skip_datasets]

    cohort = Cohort()
    for dataset_name in dataset_names:
        dataset = cohort.create_dataset(dataset_name)
        logging.info(f'Getting samples for dataset {dataset_name}')
        sample_entries = get_metamist().get_sample_entries(dataset_name)
        for entry in sample_entries:
            dataset.add_sample(
                id=str(entry['id']),
                external_id=str(entry['external_id']),
                meta=entry.get('meta', {}),
            )

    if not cohort.get_datasets():
        msg = 'No datasets populated'
        if 'skip_samples' in get_config()['workflow']:
            msg += ' (after skipping samples)'
        if 'only_samples' in get_config()['workflow']:
            msg += ' (after picking samples)'
        raise MetamistError(msg)

    if sequencing_type := get_config()['workflow'].get('sequencing_type'):
        _populate_alignment_inputs(cohort, sequencing_type)
        _filter_sequencing_type(cohort, sequencing_type)
    _populate_analysis(cohort)
    _populate_participants(cohort)
    if get_config()['workflow'].get('read_pedigree', True):
        _populate_pedigree(cohort)
    assert cohort.get_samples()
    return cohort


def _filter_sequencing_type(cohort: Cohort, sequencing_type: str):
    """
    Filtering to the samples with only requested sequencing types.
    """
    for s in cohort.get_samples():
        if not s.seq_by_type:
            s.active = False
            continue

        if s.alignment_input_by_seq_type:
            avail_types = list(s.seq_by_type.keys())
            s.alignment_input_by_seq_type = {
                k: v
                for k, v in s.alignment_input_by_seq_type.items()
                if k == sequencing_type
            }
            if not bool(s.alignment_input_by_seq_type):
                logging.warning(
                    f'{s}: skipping because no inputs with data type '
                    f'"{sequencing_type}" found in {avail_types}'
                )
                s.active = False


def _populate_alignment_inputs(
    cohort: Cohort,
    sequencing_type: str,
    check_existence: bool = False,
) -> None:
    """
    Populate sequencing inputs for samples.
    """
    assert cohort.get_sample_ids()

    seq_entries_by_sid = get_metamist().get_sequence_entries_by_sid(
        cohort.get_sample_ids(), sequencing_type=sequencing_type
    )

    # Log sequences without samples, this is a pretty common thing,
    # but useful to log to easier track down samples not processed
    if sids_wo_seq := [
        sid for sid in cohort.get_sample_ids() if sid not in seq_entries_by_sid
    ]:
        logging.info(
            f'No {sequencing_type} sequencing data found for '
            f'{len(sids_wo_seq)}/{len(cohort.get_samples())} samples:'
        )
        for ds in cohort.get_datasets():
            ds_sids_wo_seq = [sid for sid in sids_wo_seq if sid in ds.get_sample_ids()]
            logging.info(
                f'\t{ds.name}, {len(ds_sids_wo_seq)}/{len(ds.get_samples())} samples: '
                f'{", ".join(ds_sids_wo_seq)}'
            )

    sid_wo_reads = set()
    for sample in cohort.get_samples():
        for entry in seq_entries_by_sid.get(sample.id, []):
            seq = Sequence.parse(entry, parse_reads=False)
            if seq.sequencing_type in sample.seq_by_type:
                raise MetamistError(
                    f'{sample}: found more than one associated sequencing entry with '
                    f'sequencing type: {seq.sequencing_type}. Make sure there is only '
                    f'one data source of sequencing type per sample.'
                )
            sample.seq_by_type[seq.sequencing_type] = seq

            if not entry.get('meta', {}).get('reads'):
                sid_wo_reads.add(sample.id)
                continue

            alignment_input = Sequence.parse_reads(
                sample_id=sample.id,
                meta=entry['meta'],
                check_existence=check_existence,
            )
            seq.alignment_input = alignment_input
            sample.alignment_input_by_seq_type[seq.sequencing_type] = alignment_input

    if sid_wo_reads:
        logging.warning(
            f'Found {len(sid_wo_reads)}/{len(cohort.get_samples())} samples with '
            f'no meta/reads in corresponding sequence entries'
        )
        for ds in cohort.get_datasets():
            ds_sid_wo_reads = [
                sid for sid in sid_wo_reads if sid in ds.get_sample_ids()
            ]
            logging.warning(
                f'\t{ds.name}, {len(ds_sid_wo_reads)}/{len(ds.get_samples())} samples: '
                f'{", ".join(ds_sid_wo_reads)}'
            )


def _populate_analysis(cohort: Cohort) -> None:
    """
    Populate Analysis entries.
    """
    for dataset in cohort.get_datasets():
        gvcf_by_sid = get_metamist().get_analyses_by_sid(
            dataset.get_sample_ids(),
            analysis_type=AnalysisType.GVCF,
            dataset=dataset.name,
        )
        cram_by_sid = get_metamist().get_analyses_by_sid(
            dataset.get_sample_ids(),
            analysis_type=AnalysisType.CRAM,
            dataset=dataset.name,
        )
        for sample in dataset.get_samples():
            if (analysis := gvcf_by_sid.get(sample.id)) and analysis.output:
                assert analysis.output == sample.make_gvcf_path().path, (
                    analysis.output,
                    sample.make_gvcf_path().path,
                )
                sample.gvcf = sample.make_gvcf_path()
            elif sample.make_gvcf_path().exists():
                sample.gvcf = sample.make_gvcf_path()
            if (analysis := cram_by_sid.get(sample.id)) and analysis.output:
                assert analysis.output == sample.make_cram_path().path, (
                    analysis.output,
                    sample.make_cram_path().path,
                )
                sample.cram = sample.make_cram_path()
            elif sample.make_cram_path().exists():
                sample.cram = sample.make_cram_path()


def _populate_participants(cohort: Cohort) -> None:
    """
    Populate Participant entries.
    """
    for dataset in cohort.get_datasets():
        logging.info(f'Reading participants IDs for dataset {dataset}')

        participant_entry_by_sid = get_metamist().get_participant_entries_by_sid(
            dataset.name
        )

        for sample in dataset.get_samples():
            if entry := participant_entry_by_sid.get(sample.id):
                sample.participant_id = entry['external_id']
                if reported_sex := entry.get('reported_sex'):
                    sample.pedigree.sex = Sex.parse(reported_sex)
                update_dict(sample.meta, entry.get('meta', {}))


def _populate_pedigree(cohort: Cohort) -> None:
    """
    Populate pedigree data for samples.
    """
    sample_by_participant_id = dict()
    for s in cohort.get_samples():
        sample_by_participant_id[s.participant_id] = s

    for dataset in cohort.get_datasets():
        logging.info(f'Reading pedigree for dataset {dataset}')
        ped_entries = get_metamist().get_ped_entries(dataset=dataset.name)
        ped_entry_by_participant_id = {}
        for ped_entry in ped_entries:
            part_id = str(ped_entry['individual_id'])
            ped_entry_by_participant_id[part_id] = ped_entry

        sids_wo_ped = []
        for sample in dataset.get_samples():
            if sample.participant_id not in ped_entry_by_participant_id:
                sids_wo_ped.append(sample.id)
                continue

            ped_entry = ped_entry_by_participant_id[sample.participant_id]
            maternal_sample = sample_by_participant_id.get(
                str(ped_entry['maternal_id'])
            )
            paternal_sample = sample_by_participant_id.get(
                str(ped_entry['paternal_id'])
            )
            sample.pedigree = PedigreeInfo(
                sample=sample,
                fam_id=ped_entry['family_id'],
                mom=maternal_sample,
                dad=paternal_sample,
                sex=Sex.parse(str(ped_entry['sex'])),
                phenotype=ped_entry['affected'] or '0',
            )
        if sids_wo_ped:
            logging.warning(
                f'No pedigree data found for '
                f'{len(sids_wo_ped)}/{len(dataset.get_samples())} samples'
            )
