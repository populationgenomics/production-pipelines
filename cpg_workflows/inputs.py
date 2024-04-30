"""
Metamist wrapper to get input sequencing groups.
"""

import logging

from cpg_utils.config import get_config, update_dict
from cpg_workflows.filetypes import CramPath, GvcfPath

from .metamist import AnalysisType, Assay, MetamistError, get_metamist, parse_reads
from .targets import Cohort, PedigreeInfo, SequencingGroup, Sex
from .utils import exists

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
    dataset_names = get_config()['workflow'].get('input_datasets', [analysis_dataset_name])
    skip_datasets = get_config()['workflow'].get('skip_datasets', [])
    dataset_names = [d for d in dataset_names if d not in skip_datasets]

    cohort = Cohort()
    for dataset_name in dataset_names:
        dataset = cohort.create_dataset(dataset_name)
        sequencing_group_entries = get_metamist().get_sg_entries(dataset_name)
        for entry in sequencing_group_entries:
            metadata = entry.get('meta', {})
            update_dict(metadata, entry['sample']['participant'].get('meta', {}))

            # phenotypes are managed badly here, need a cleaner way to get them into the SG
            update_dict(metadata, {'phenotypes': entry['sample']['participant'].get('phenotypes', {})})

            sequencing_group = dataset.add_sequencing_group(
                id=str(entry['id']),
                external_id=str(entry['sample']['externalId']),
                participant_id=entry['sample']['participant'].get('externalId'),
                meta=metadata,
            )

            if reported_sex := entry['sample']['participant'].get('reportedSex'):
                sequencing_group.pedigree.sex = Sex.parse(reported_sex)

            _populate_alignment_inputs(sequencing_group, entry)

    if not cohort.get_datasets():
        msg = 'No datasets populated'
        if 'skip_sgs' in get_config()['workflow']:
            msg += ' (after skipping sequencing groups)'
        if 'only_sgs' in get_config()['workflow']:
            msg += ' (after picking sequencing groups)'
        raise MetamistError(msg)

    _populate_analysis(cohort)
    if get_config()['workflow'].get('read_pedigree', True):
        _populate_pedigree(cohort)
    assert cohort.get_sequencing_groups()
    return cohort


def _combine_assay_meta(assays: list[Assay]) -> dict:
    """
    Combine assay meta from multiple assays
    """
    assays_meta: dict = {}
    for assay in assays:
        for key, value in assay.meta.items():
            if key == 'reads':
                assays_meta.setdefault(key, []).append(value)
            else:
                assays_meta[key] = value
    return assays_meta


def _populate_alignment_inputs(
    sequencing_group: SequencingGroup,
    entry: dict,
    check_existence: bool = False,
) -> None:
    """
    Populate sequencing inputs for a sequencing group
    """
    assays: list[Assay] = []
    for assay in entry['assays']:
        _assay = Assay.parse(assay, sequencing_group.id, run_parse_reads=False)
        assays.append(_assay)

    # Check only one assay type per sequencing group
    assert len(set([assay.assay_type for assay in assays])) == 1
    sequencing_group.assays[assays[0].assay_type] = tuple(assays)

    if len(entry['assays']) > 1:
        _assay_meta = _combine_assay_meta(assays)
    else:
        _assay_meta = assays[0].meta
        if _assay_meta.get('reads'):
            _assay_meta['reads'] = [_assay_meta['reads']]

    if _assay_meta.get('reads'):
        alignment_input = parse_reads(
            sequencing_group_id=sequencing_group.id,
            assay_meta=_assay_meta,
            check_existence=check_existence,
        )
        sequencing_group.alignment_input_by_seq_type[entry['type']] = alignment_input
    else:
        logging.warning(f'No reads found for sequencing group {sequencing_group.id} of type {entry["type"]}')

    return None


def _populate_analysis(cohort: Cohort) -> None:
    """
    Populate Analysis entries.
    """
    for dataset in cohort.get_datasets():
        gvcf_by_sgid = get_metamist().get_analyses_by_sgid(
            dataset.get_sequencing_group_ids(),
            analysis_type=AnalysisType.GVCF,
            dataset=dataset.name,
        )
        cram_by_sgid = get_metamist().get_analyses_by_sgid(
            dataset.get_sequencing_group_ids(),
            analysis_type=AnalysisType.CRAM,
            dataset=dataset.name,
        )

        for sequencing_group in dataset.get_sequencing_groups():
            if (analysis := gvcf_by_sgid.get(sequencing_group.id)) and analysis.output:
                # assert file exists
                assert exists(analysis.output), (
                    'gvcf file does not exist',
                    analysis.output,
                )
                sequencing_group.gvcf = GvcfPath(path=analysis.output)
            elif exists(sequencing_group.make_gvcf_path()):
                logging.warning(
                    f'We found a gvcf file in the expected location {sequencing_group.make_gvcf_path()},'
                    'but it is not logged in metamist. Skipping. You may want to update the metadata and try again. ',
                )
            if (analysis := cram_by_sgid.get(sequencing_group.id)) and analysis.output:
                # assert file exists
                assert exists(analysis.output), (
                    'cram file does not exist',
                    analysis.output,
                )
                crai_path = analysis.output.with_suffix('.cram.crai')
                if not exists(crai_path):
                    crai_path = None
                sequencing_group.cram = CramPath(analysis.output, crai_path)

            elif exists(sequencing_group.make_cram_path()):
                logging.warning(
                    f'We found a cram file in the expected location {sequencing_group.make_cram_path()},'
                    'but it is not logged in metamist. Skipping. You may want to update the metadata and try again. ',
                )


def _populate_pedigree(cohort: Cohort) -> None:
    """
    Populate pedigree data for sequencing groups.
    """
    sg_by_participant_id = dict()
    for sg in cohort.get_sequencing_groups():
        sg_by_participant_id[sg.participant_id] = sg

    for dataset in cohort.get_datasets():
        logging.info(f'Reading pedigree for dataset {dataset}')
        ped_entries = get_metamist().get_ped_entries(dataset=dataset.name)
        ped_entry_by_participant_id = {}
        for ped_entry in ped_entries:
            part_id = str(ped_entry['individual_id'])
            ped_entry_by_participant_id[part_id] = ped_entry

        sgids_wo_ped = []
        for sequencing_group in dataset.get_sequencing_groups():
            if sequencing_group.participant_id not in ped_entry_by_participant_id:
                sgids_wo_ped.append(sequencing_group.id)
                continue

            ped_entry = ped_entry_by_participant_id[sequencing_group.participant_id]
            maternal_sg = sg_by_participant_id.get(str(ped_entry['maternal_id']))
            paternal_sg = sg_by_participant_id.get(str(ped_entry['paternal_id']))
            sequencing_group.pedigree = PedigreeInfo(
                sequencing_group=sequencing_group,
                fam_id=ped_entry['family_id'],
                mom=maternal_sg,
                dad=paternal_sg,
                sex=Sex.parse(str(ped_entry['sex'])),
                phenotype=ped_entry['affected'] or '0',
            )
        if sgids_wo_ped:
            logging.warning(
                f'No pedigree data found for '
                f'{len(sgids_wo_ped)}/{len(dataset.get_sequencing_groups())} sequencing groups',
            )
