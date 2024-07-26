"""
Metamist wrapper to get input sequencing groups.
"""

import logging

from cpg_utils.config import config_retrieve, update_dict
from cpg_workflows.filetypes import CramPath, GvcfPath

from .metamist import AnalysisType, Assay, MetamistError, get_metamist, parse_reads
from .targets import Cohort, MultiCohort, PedigreeInfo, SequencingGroup, Sex
from .utils import exists

_cohort: Cohort | None = None
_multicohort: MultiCohort | None = None


def actual_get_multicohort() -> MultiCohort:
    """Return the multicohort object"""
    global _multicohort
    if not _multicohort:
        _multicohort = create_multicohort()
    return _multicohort


def deprecated_get_cohort() -> Cohort:
    """
    This should still return a MultiCohort, but with a single Cohort
    i.e. in this case Cohort and MultiCohort are the same
    """
    global _multicohort
    if not _multicohort:
        _multicohort = deprecated_create_multicohort()
    return _multicohort


def get_multicohort() -> Cohort | MultiCohort:
    """
    Return the cohort or multicohort object based on the workflow configuration.
    """
    input_datasets = config_retrieve(['workflow', 'input_datasets'], [])
    custom_cohort_ids = config_retrieve(['workflow', 'input_cohorts'], [])
    if custom_cohort_ids and input_datasets:
        raise ValueError('Cannot use both custom_cohort_ids and input_datasets in the same workflow')

    # NOTE: When configuring sgs in the config is deprecated, this will be removed.
    if config_retrieve(['workflow', 'input_cohorts'], []):
        if not isinstance(custom_cohort_ids, list):
            raise ValueError('input_cohorts must be a list')
        return actual_get_multicohort()
    return deprecated_get_cohort()


def create_multicohort() -> MultiCohort:
    """
    Add cohorts in the multicohort.
    """
    config = config_retrieve(['workflow'])
    multicohort = MultiCohort()
    read_pedigree = config.get('read_pedigree', True)

    for cohort_id in config.get('input_cohorts', []):
        # create a shell for this Cohort
        cohort = multicohort.create_or_return_cohort(cohort_id)

        # get a list of sequencing group data dictionaries for the cohort
        sgs_in_cohort = get_metamist().get_sgs_from_cohort(cohort_id)

        # populate the Cohort with SequencingGroups
        _populate_cohort(cohort, sgs=sgs_in_cohort, read_pedigree=read_pedigree)

    for sg in multicohort.get_sequencing_groups():
        multicohort.create_or_return_dataset(sg.dataset).add_sequencing_group(sg)

    # now populate the per-dataset analyses
    for dataset in multicohort.get_datasets():
        _populate_analysis(dataset)
        if read_pedigree:
            _populate_pedigree(dataset)

    return multicohort


def _populate_cohort(cohort: Cohort, sgs: list[dict], read_pedigree: bool = True):
    """
    Add sequencing groups to the cohort
    Multiple cohorts can co-exist in the workflow run, and could have overlapping SG ID sets
    """

    for entry in sgs:
        metadata = entry.get('meta', {})
        update_dict(metadata, entry['sample']['participant'].get('meta', {}))
        # phenotypes are managed badly here, need a cleaner way to get them into the SG
        update_dict(metadata, {'phenotypes': entry['sample']['participant'].get('phenotypes', {})})

        sequencing_group = SequencingGroup(
            id=entry['id'],
            external_id=str(entry['sample']['externalId']),
            dataset=entry['sample']['project']['name'],
            participant_id=entry['sample']['participant'].get('externalId'),
            meta=metadata,
            sequencing_type=entry['type'],
            sequencing_technology=entry['technology'],
            sequencing_platform=entry['platform'],
        )

        if reported_sex := entry['sample']['participant'].get('reportedSex'):
            sequencing_group.pedigree.sex = Sex.parse(reported_sex)

        _populate_alignment_inputs(sequencing_group, entry)

        # add this SequencingGroup to the Cohort
        cohort.add_sequencing_group(sequencing_group)

    assert cohort.get_sequencing_groups()


def deprecated_create_multicohort() -> MultiCohort:
    """
    Create a MultiCohort from component input datasets
    The MultiCohort will contain one or more Datasets
    and a single Cohort containing all SequencingGroups
    """

    # create the MultiCohort shell
    multi_cohort = MultiCohort()

    # One cohort for the whole workflow run (equal to the MultiCohort)
    cohort = Cohort()

    config = config_retrieve(['workflow'])
    analysis_dataset_name = config_retrieve(['workflow', 'dataset'])
    input_datasets = config_retrieve(['workflow', 'input_datasets'], [])
    skip_datasets = config_retrieve(['workflow', 'skip_datasets'], [])

    if input_datasets:
        dataset_names = input_datasets
        logging.warning('Using input_datasets will soon be deprecated. Use input_cohorts instead.')
    else:
        dataset_names = [analysis_dataset_name]
        logging.warning('Using dataset will soon be deprecated. Use input_cohorts instead.')

    dataset_names = [d for d in dataset_names if d not in skip_datasets]

    for dataset_name in dataset_names:
        # dataset = cohort.create_dataset(dataset_name)
        sgs = get_metamist().get_sg_entries(dataset_name)
        dataset = Dataset(name=dataset_name, multicohort=multi_cohort)

        for entry in sgs:
            metadata = entry.get('meta', {})
            update_dict(metadata, entry['sample']['participant'].get('meta', {}))
            # phenotypes are managed badly here, need a cleaner way to get them into the SG
            update_dict(metadata, {'phenotypes': entry['sample']['participant'].get('phenotypes', {})})

            # create the SequencingGroup object
            sequencing_group = SequencingGroup(
                id=entry['id'],
                external_id=str(entry['sample']['externalId']),
                dataset=dataset_name,
                participant_id=entry['sample']['participant'].get('externalId'),
                meta=metadata,
                sequencing_type=entry['type'],
                sequencing_technology=entry['technology'],
                sequencing_platform=entry['platform'],
            )

            if reported_sex := entry['sample']['participant'].get('reportedSex'):
                sequencing_group.pedigree.sex = Sex.parse(reported_sex)

            _populate_alignment_inputs(sequencing_group, entry)

            # add this SequencingGroup to the Dataset
            dataset.add_sequencing_group(sequencing_group)
            # also add this SequencingGroup to the Cohort
            cohort.add_sequencing_group(sequencing_group)

        # populate the analyses for this dataset
        _populate_analysis(dataset)
        multi_cohort.add_dataset_to_multicohort(dataset)

    if not cohort.get_datasets():
        msg = 'No datasets populated'
        if 'skip_sgs' in config:
            msg += ' (after skipping sequencing groups)'
        if 'only_sgs' in config:
            msg += ' (after picking sequencing groups)'
        raise MetamistError(msg)

    if config.get('read_pedigree', True):
        _populate_pedigree(cohort)

    multi_cohort.add_cohort_to_multicohort(cohort)

    return multi_cohort


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
    should this be a method of SequencingGroup?
    """
    assays: list[Assay] = []

    for assay in entry.get('assays', []):
        _assay = Assay.parse(assay, sequencing_group.id, run_parse_reads=False)
        assays.append(_assay)

    sequencing_group.assays = tuple(assays)

    if len(assays) > 1:
        _assay_meta = _combine_assay_meta(assays)
    else:
        _assay_meta = assays[0].meta
        if _reads := _assay_meta.get('reads'):
            _assay_meta['reads'] = [_reads]

    if _assay_meta.get('reads'):
        alignment_input = parse_reads(
            sequencing_group_id=sequencing_group.id,
            assay_meta=_assay_meta,
            check_existence=check_existence,
        )
        sequencing_group.alignment_input = alignment_input
    else:
        logging.warning(f'No reads found for sequencing group {sequencing_group.id} of type {entry["type"]}')


def _populate_analysis(dataset: Dataset) -> None:
    """
    Populate Analysis entries.
    # todo: we should postpone calling this method until we gathered all SG IDs in all cohorts,
    # todo: then binned the resulting SG IDs by dataset - call this method once per dataset
    """
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


def _populate_pedigree(dataset: Dataset) -> None:
    """
    Populate pedigree data for sequencing groups.
    """

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
        maternal_sg = dataset.get_sequencing_group_by_id(str(ped_entry['maternal_id']))
        paternal_sg = dataset.get_sequencing_group_by_id(str(ped_entry['paternal_id']))
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
