"""
Metamist wrapper to get input samples.
"""

import logging

from cpg_utils.config import get_config, update_dict

from .metamist import get_metamist, Assay, AnalysisType, MetamistError
from .targets import Cohort, Sex, PedigreeInfo, SequencingGroup


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
        logging.info(f'Getting sequencing groups for dataset {dataset_name}')
        sequencing_group_entries = get_metamist().get_sg_entries(dataset_name)
        for entry in sequencing_group_entries:
            # TODO: Add participant external ID here
            # Add pedigree sex.
            metadata = entry.get('meta', {})
            update_dict(metadata, entry['sample']['participant'].get('meta', {}))

            sequencing_group = dataset.add_sequencing_group(
                id=str(entry['id']),
                external_id=str(entry['sample']['externalId']),
                participant_id=entry['sample']['participant'].get('externalId'),
                meta=metadata,
            )

            if reported_sex := entry['sample']['participant'].get('reportedSex'):
                sequencing_group.pedigree.sex = Sex.parse(reported_sex)

            _populate_alignment_inputs(sequencing_group, entry)
            # TODO: Add checks here,logging samples without sequences
            # Check there is only one sequencing group per type.

    if not cohort.get_datasets():
        msg = 'No datasets populated'
        if 'skip_sgs' in get_config()['workflow']:
            msg += ' (after skipping sequencing groups)'
        if 'only_sgs' in get_config()['workflow']:
            msg += ' (after picking sequencing groups)'
        raise MetamistError(msg)

    # NOTE: Are there cases where there isn't a sequencing_type?
    _populate_analysis(cohort)
    if get_config()['workflow'].get('read_pedigree', True):
        _populate_pedigree(cohort)
    assert cohort.get_sequencing_groups()
    return cohort


def _populate_alignment_inputs(
    sequencing_group: SequencingGroup,
    entry: dict,
    check_existence: bool = False,
) -> None:
    """
    Populate sequencing inputs for a sequencing group
    """

    if len(entry['assays']) != 1:
        raise MetamistError(
            f'Only one active assay per sequencing group allowed. Found {len(entry["assays"])}'
        )

    assay_entry = entry['assays'][0]

    assay = Assay.parse(assay_entry, sequencing_group.id, parse_reads=False)

    sequencing_group.assays[assay.assay_type] = assay

    if assay_entry.get('meta', {}).get('reads'):
        alignment_input = Assay.parse_reads(
            sample_id=sequencing_group.id,
            meta=assay_entry['meta'],
            check_existence=check_existence,
        )

    assay.alignment_input = alignment_input
    sequencing_group.alignment_input_by_seq_type[entry['type']] = alignment_input

    # TODO: Include some additional logging here about sequences without reads etc.

    return None


def _populate_analysis(cohort: Cohort) -> None:
    """
    Populate Analysis entries.
    """
    for dataset in cohort.get_datasets():
        gvcf_by_sid = get_metamist().get_analyses_by_sid(
            dataset.get_sequencing_group_ids(),
            analysis_type=AnalysisType.GVCF,
            dataset=dataset.name,
        )
        cram_by_sid = get_metamist().get_analyses_by_sid(
            dataset.get_sequencing_group_ids(),
            analysis_type=AnalysisType.CRAM,
            dataset=dataset.name,
        )

        # NOTE: This logic will be simplified to remove existence checks overwriting metamist
        # in a later PR.
        for sequencing_group in dataset.get_sequencing_groups():
            if (analysis := gvcf_by_sid.get(sequencing_group.id)) and analysis.output:
                assert analysis.output == sequencing_group.make_gvcf_path().path, (
                    analysis.output,
                    sequencing_group.make_gvcf_path().path,
                )
                sequencing_group.gvcf = sequencing_group.make_gvcf_path()
            elif sequencing_group.make_gvcf_path().exists():
                sequencing_group.gvcf = sequencing_group.make_gvcf_path()
            if (analysis := cram_by_sid.get(sequencing_group.id)) and analysis.output:
                assert analysis.output == sequencing_group.make_cram_path().path, (
                    analysis.output,
                    sequencing_group.make_cram_path().path,
                )
                sequencing_group.cram = sequencing_group.make_cram_path()
            elif sequencing_group.make_cram_path().exists():
                sequencing_group.cram = sequencing_group.make_cram_path()


def _populate_pedigree(cohort: Cohort) -> None:
    """
    Populate pedigree data for samples.
    """
    sample_by_participant_id = dict()
    for s in cohort.get_sequencing_groups():
        sample_by_participant_id[s.participant_id] = s

    for dataset in cohort.get_datasets():
        logging.info(f'Reading pedigree for dataset {dataset}')
        ped_entries = get_metamist().get_ped_entries(dataset=dataset.name)
        ped_entry_by_participant_id = {}
        for ped_entry in ped_entries:
            part_id = str(ped_entry['individual_id'])
            ped_entry_by_participant_id[part_id] = ped_entry

        sids_wo_ped = []
        for sequencing_group in dataset.get_sequencing_groups():
            if sequencing_group.participant_id not in ped_entry_by_participant_id:
                sids_wo_ped.append(sequencing_group.id)
                continue

            ped_entry = ped_entry_by_participant_id[sequencing_group.participant_id]
            maternal_sample = sample_by_participant_id.get(
                str(ped_entry['maternal_id'])
            )
            paternal_sample = sample_by_participant_id.get(
                str(ped_entry['paternal_id'])
            )
            sequencing_group.pedigree = PedigreeInfo(
                sequencing_group=sequencing_group,
                fam_id=ped_entry['family_id'],
                mom=maternal_sample,
                dad=paternal_sample,
                sex=Sex.parse(str(ped_entry['sex'])),
                phenotype=ped_entry['affected'] or '0',
            )
        if sids_wo_ped:
            logging.warning(
                f'No pedigree data found for '
                f'{len(sids_wo_ped)}/{len(dataset.get_sequencing_groups())} sequencing groups'
            )
