"""
Metamist wrapper to get input samples.
"""

import logging
import traceback
from collections import defaultdict
from itertools import groupby

from cpg_utils.config import get_config
from sample_metadata import ApiException

from cpg_pipes.metamist import get_metamist, MmSequence
from cpg_pipes.targets import Cohort, Sex, PedigreeInfo

logger = logging.getLogger(__file__)


class InputError(Exception):
    """
    Exception thrown when there is something wrong while parsing inputs.
    """

    pass


_cohort: Cohort | None = None


def get_cohort() -> Cohort:
    """Return the cohort object, whic is a signleton"""
    global _cohort
    if not _cohort:
        _cohort = _create_cohort()
    return _cohort


def _create_cohort() -> Cohort:
    """
    Add datasets in the cohort. There exists only one cohort for
    the pipeline run.
    """
    analysis_dataset_name = get_config()['workflow']['dataset']
    dataset_names = get_config()['workflow'].get('datasets', [analysis_dataset_name])
    skip_datasets = get_config()['workflow'].get('skip_datasets', [])
    dataset_names = [d for d in dataset_names if d not in skip_datasets]

    skip_samples = get_config()['workflow'].get('skip_samples', [])
    only_samples = get_config()['workflow'].get('only_samples', [])

    cohort = Cohort()
    for dataset_name in dataset_names:
        dataset = cohort.create_dataset(dataset_name)
        sample_entries = get_metamist().sapi.get_samples(
            body_get_samples={'project_ids': [dataset_name]}
        )
        sample_entries = _filter_samples(
            sample_entries,
            dataset_name,
            skip_samples,
            only_samples,
        )
        for entry in sample_entries:
            dataset.add_sample(
                id=str(entry['id']),
                external_id=str(entry['external_id']),
                meta=entry['meta'],
            )

    if not cohort.get_datasets():
        msg = 'No datasets populated'
        if skip_samples or only_samples or skip_datasets:
            msg += ' (after skipping/picking samples)'
        logger.warning(msg)
        return cohort

    _populate_alignment_inputs(cohort)
    _filter_sequencing_type(cohort)
    _populate_analysis(cohort)
    _populate_participants(cohort)
    _populate_pedigree(cohort)
    return cohort


def _filter_sequencing_type(cohort: Cohort):
    """
    Filtering to the samples with only requested sequencing types.
    """
    sequencing_type = get_config()['workflow']['sequencing_type']
    for s in cohort.get_samples():
        if not s.alignment_input_by_seq_type:
            logger.warning(f'{s}: skipping because no sequencing inputs found')
            s.active = False
            continue

        avail_types = list(s.alignment_input_by_seq_type.keys())
        s.alignment_input_by_seq_type = {
            k: v
            for k, v in s.alignment_input_by_seq_type.items()
            if k == sequencing_type
        }
        if not bool(s.alignment_input_by_seq_type):
            logger.warning(
                f'{s}: skipping because no inputs with data type '
                f'"{sequencing_type}" found in {avail_types}'
            )
            s.active = False


def _filter_samples(
    entries: list[dict[str, str]],
    dataset_name: str,
    skip_samples: list[str] | None = None,
    only_samples: list[str] | None = None,
) -> list[dict]:
    """
    Apply the only_samples and skip_samples filters.
    """

    filtered_entries = []
    for entry in entries:
        cpgid = entry['id']
        extid = entry['external_id']
        if only_samples:
            if cpgid in only_samples or extid in only_samples:
                logger.info(f'Picking sample: {dataset_name}|{cpgid}|{extid}')
            else:
                continue
        if skip_samples:
            if cpgid in skip_samples or extid in skip_samples:
                logger.info(f'Skipping sample: {dataset_name}|{cpgid}|{extid}')
                continue
        filtered_entries.append(entry)
    return filtered_entries


def _populate_alignment_inputs(
    cohort: Cohort,
    check_existence: bool = False,
) -> None:
    """
    Populate sequencing inputs for samples.
    """
    assert cohort.get_sample_ids()
    seq_type = get_config()['workflow']['sequencing_type']
    try:
        found_seqs: list[dict] = get_metamist().seqapi.get_sequences_by_sample_ids(
            cohort.get_sample_ids(), get_latest_sequence_only=False
        )
        found_seqs = [seq for seq in found_seqs if str(seq['type']) == seq_type]
    except ApiException:
        if get_config()['workflow'].get('smdb_errors_are_fatal', True):
            raise
        else:
            logger.error(
                'Getting sequencing data from SMDB resulted in an error. '
                'Continuing without sequencing data because of flag override. '
                'However, here is the error: '
            )
            traceback.print_exc()
            found_seqs = []
    found_seqs_by_sid = defaultdict(list)
    for found_seq in found_seqs:
        found_seqs_by_sid[found_seq['sample_id']].append(found_seq)

    # Log sequences without samples, this is a pretty common thing,
    # but useful to log to easier track down samples not processed
    if sample_wo_seq := [
        s for s in cohort.get_samples() if s.id not in found_seqs_by_sid
    ]:
        msg = f'No {seq_type} sequencing data found for samples:\n'
        ds_sample_count = {
            ds_name: len(list(ds_samples))
            for ds_name, ds_samples in groupby(
                cohort.get_samples(), key=lambda s: s.dataset.name
            )
        }
        for ds, samples in groupby(sample_wo_seq, key=lambda s: s.dataset.name):
            msg += (
                f'\t{ds}, {len(list(samples))}/{ds_sample_count.get(ds)} samples: '
                f'{", ".join([s.id for s in samples])}\n'
            )
        logger.info(msg)

    for sample in cohort.get_samples():
        for d in found_seqs_by_sid.get(sample.id, []):
            seq = MmSequence.parse(d, check_existence=check_existence)
            if seq.alignment_input:
                if seq.sequencing_type in sample.alignment_input_by_seq_type:
                    raise InputError(
                        f'{sample}: found more than 1 alignment input with '
                        f'sequencing type: {seq.sequencing_type}. Check your '
                        f'input provider to make sure there is only one data source '
                        f'of sequencing type per sample.'
                    )
                sample.alignment_input_by_seq_type[
                    seq.sequencing_type
                ] = seq.alignment_input


def _populate_analysis(cohort: Cohort) -> None:
    """
    Populate Analysis entries.
    """
    pass


def _populate_participants(cohort: Cohort) -> None:
    """
    Populate Participant entries.
    """
    for dataset in cohort.get_datasets():
        pid_sid_multi = (
            get_metamist().papi.get_external_participant_id_to_internal_sample_id(
                dataset.name
            )
        )
        participant_by_sid = {}
        for group in pid_sid_multi:
            pid = group[0]
            for sid in group[1:]:
                participant_by_sid[sid] = pid.strip()

        for sample in dataset.get_samples():
            if pid := participant_by_sid.get(sample.id):
                sample.participant_id = pid


def _populate_pedigree(cohort: Cohort) -> None:
    """
    Populate pedigree data for samples.
    """
    sample_by_participant_id = dict()
    for s in cohort.get_samples():
        sample_by_participant_id[s.participant_id] = s

    for dataset in cohort.get_datasets():
        ped_entries = get_metamist().get_ped_entries(dataset=dataset.name)
        ped_entry_by_participant_id = {}
        for ped_entry in ped_entries:
            part_id = str(ped_entry['individual_id'])
            ped_entry_by_participant_id[part_id] = ped_entry

        for sample in dataset.get_samples():
            if sample.participant_id not in ped_entry_by_participant_id:
                logger.warning(
                    f'No pedigree data for participant {sample.participant_id}'
                )
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

    for dataset in cohort.get_datasets():
        samples_with_ped = [s for s in dataset.get_samples() if s.pedigree]
        logger.info(
            f'{dataset.name}: found pedigree info for {len(samples_with_ped)} '
            f'samples out of {len(dataset.get_samples())}'
        )
