#!/usr/bin/env python3
# pylint: disable=too-many-instance-attributes,too-many-lines,too-many-locals

""" Example Invocation

analysis-runner \
--dataset acute-care --description "populate acute care test subset" --output-dir "acute-care-test" \
--access-level full \
scripts/create_test_subset.py --project acute-care --families 4

This example will populate acute-care-test with the metamist data for 4 families.
"""

import csv
import logging
import os
import random
import subprocess
from argparse import ArgumentParser, BooleanOptionalAction
from collections import Counter
from typing import Any, Tuple

from google.cloud import storage

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

from metamist.apis import AnalysisApi, FamilyApi, ParticipantApi, SampleApi
from metamist.graphql import gql, query
from metamist.models import (
    Analysis,
    AnalysisStatus,
    AnalysisUpdateModel,
    AssayUpsert,
    SampleUpsert,
    SequencingGroupUpsert,
)

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

sapi = SampleApi()
aapi = AnalysisApi()
fapi = FamilyApi()
papi = ParticipantApi()

DEFAULT_SAMPLES_N = 10

PRIMARY_EXTERNAL_ORG = ''

QUERY_ALL_DATA = gql(
    """
    query getAllData($project: String!, $sids: [String!]) {
        project(name: $project) {
            samples(id: {in_: $sids}) {
                id
                meta
                type
                externalId
                externalIds
                participant {
                    externalId
                    externalIds
                    id
                    karyotype
                    meta
                    reportedGender
                    reportedSex
                }
                sequencingGroups{
                    id
                    meta
                    platform
                    technology
                    type
                    assays {
                        id
                        meta
                        type
                    }
                    analyses(type: {in_: ["cram", "gvcf"]}) {
                        active
                        id
                        meta
                        output
                        status
                        timestampCompleted
                        type
                    }
                }
            }
        }
    }
    """
)

# TODO: We can change this to filter external sample ids
EXISTING_DATA_QUERY = gql(
    """
    query getExistingData($project: String!) {
        project(name: $project) {
            samples{
                id
                externalId
                sequencingGroups {
                    id
                    type
                    assays {
                        id
                        meta
                        type
                    }
                    analyses {
                        id
                        type
                    }
                }
            }
        }
    }
    """
)

QUERY_FAMILY_SAMPLES = gql(
    """
    query FamilyQuery($project: String!) {
        project(name: $project) {
            families {
                id
                externalId
                participants {
                    samples {
                        id
                    }
                }
            }
        }
    }
"""
)

SG_ID_QUERY = gql(
    """
    query getSGIds($project: String!) {
        project(name: $project) {
            samples{
                id
                externalId
                sequencingGroups {
                    id
                    type
                    platform
                    technology
                }
            }
        }
    }
    """
)

PARTICIPANT_QUERY = gql(
    """
    query ($project: String!) {
        project (name: $project) {
            participants {
                externalIds
                id
                meta
                phenotypes
                reportedGender
                reportedSex
            }
        }
    }
    """
)

COHORT_QUERY = gql(
    """
    query CohortQuery($project: String!) {
        project(name: $project) {
            cohorts {
                id
                sequencingGroups {
                    sample {
                        id
                    }
                }
            }
        }
    }
    """
)

# take a parameter which is an array of ints
QUERY_FAMILY_BY_INTERNAL_IDS = gql(
    """
    query FamilyQuery($project: String!, $internal_ids: [Int!]!) {
  project(name: $project) {
    participants(id: {in_: $internal_ids}) {
      families {
        id
        codedPhenotype
        description
        externalIds
      }
    }
  }
}
"""
)


def get_sids_by_random_sampling(project: str, sample_count: int, samples_so_far: set[str] | None = None) -> list[str]:
    """
    Take a project to select samples from, a number of samples to select, and a set of samples already selected.
    Select sample_count new samples from the project, excluding any already selected.

    Args:
        project (str): the metamist project to query
        sample_count (int): number of additional samples to add
        samples_so_far (set[str] | None): optional, samples we already selected by other means, remove before sampling

    Returns:
        list[str]: a list of sample IDs
    """

    if samples_so_far is None:
        samples_so_far = set()

    logger.info(f'Querying all samples in {project}')
    sid_output = query(SG_ID_QUERY, variables={'project': project})

    # all_sids here is a set of all internal IDs in the project, minus any we already selected
    all_sids = {sid['id'] for sid in sid_output.get('project').get('samples')} - samples_so_far

    logger.info(f'Found {len(all_sids)} sample ids in {project}, excluding any previously seen')

    # Randomly select from the remaining sgs
    return random.sample(sorted(all_sids), sample_count)


def main(
    project: str,
    samples_n: int,
    families_n: int,
    cohort_samples_n: int,
    additional_families: set[str],
    additional_samples: set[str],
    cohorts: set[str],
    skip_ped: bool,
    update_embedded_ids: bool,
):
    """
    Script creates a test subset for a given project.
    A new project with a suffix -test is created, and for any files in sample/meta,
    sequence/meta, or analysis/output a copy in the -test namespace is created.
    """
    if not any(
        [additional_families, additional_samples, samples_n, families_n, cohorts]
    ):
        raise ValueError('Come on, what exactly are you asking for?')

    # TODO(vbakiris) IMO it makes sense to pull the whole cohort if we're requesting a cohort?
    # TODO(vbakiris) for any smaller subset we have other params to select for that
    if cohorts and not cohort_samples_n:
        raise ValueError(
            'You must specify the number of samples to transfer from the cohort.'
        )

    # for reproducibility
    logger.info('Setting random seed to 42')
    random.seed(42)

    # 1a. Find SG IDs to be moved by Family ID to -test.
    if families_n or additional_families:
        additional_samples.update(
            get_sids_for_families(project, families_n, additional_families)
        )

    # 1b. Find SG IDs to be moved by Cohort ID to -test.
    if cohorts:
        logger.info(f'Updating sample selection to use cohorts: {cohorts}')
        additional_samples.update(
            get_sids_for_cohorts(project, cohorts, cohort_samples_n)
        )

    # 1c. If we are gathering random additional samples, get all sample IDs in project and sample from the collection
    if samples_n:
        logger.info(f'Updating sample selection with {samples_n} additional random samples')
        additional_samples.update(
            get_sids_by_random_sampling(
                project=project,
                sample_count=samples_n,
                samples_so_far=additional_samples
            )
        )

    # 2. Query all the samples from the selected sgs
    logger.info(f'Transferring {len(additional_samples)} samples. Querying metadata.')
    original_project_subset_data = query(
        QUERY_ALL_DATA, {'project': project, 'sids': list(additional_samples)}
    )

    # Pull Participant Data
    participant_data = []
    internal_participant_ids: list = []
    for sg in original_project_subset_data.get('project').get('samples', []):
        participant = sg.get('participant')
        if participant:
            participant_data.append(participant)
            internal_participant_ids.append(participant.get('id'))

    # Populating test project
    target_project = project + '-test'

    # Parse Families & Participants
    if skip_ped:
        # If no family data is available, only the participants should be transferred.
        logger.info(
            'Skipping pedigree/family information. Transferring participants only.'
        )
        logger.info(f'Transferring {len(participant_data)} participants. ')
        upserted_participant_map = transfer_participants(
            target_project=target_project,
            participant_data=participant_data,
        )

    else:
        logger.info(f'Transferring {len(participant_data)} participants. ')
        family_ids = transfer_families(
            project, target_project, internal_participant_ids
        )
        upserted_participant_map = transfer_ped(project, target_project, family_ids)

    existing_data = query(EXISTING_DATA_QUERY, {'project': target_project})

    # Create batch if needed for running header updating jobs
    if update_embedded_ids:
        this_batch = os.environ.get('HAIL_BATCH_ID')
        this_batch_text = f' invoked by batch {this_batch}' if this_batch else ''
        batch = get_batch(f'Reheadering for create_test_subset.py{this_batch_text}')
    else:
        batch = None

    logger.info('Transferring samples, sequencing groups, and assays')
    samples = original_project_subset_data.get('project').get('samples')
    old_sid_to_new_sid, sample_to_sg_attribute_map = transfer_samples_sgs_assays(
        samples=samples,
        existing_data=existing_data,
        upserted_participant_map=upserted_participant_map,
        target_project=target_project,
        project=project
    )

    logger.info('Transferring analyses')
    transfer_analyses(
        samples,
        existing_data,
        target_project,
        project,
        old_sid_to_new_sid,
        sample_to_sg_attribute_map,
        update_embedded_ids,
    )
    logger.info('Subset generation complete!')

    if batch:
        logger.info('Submitting batch to update IDs embedded in file headers')
        batch.run(wait=False)


def transfer_samples_sgs_assays(
    samples: dict,
    existing_data: dict,
    upserted_participant_map: dict[str, int],
    target_project: str,
    project: str,
):
    """
    Transfer samples, sequencing groups, and assays from the original project to the
    test project. Also creates a mapping from sample identifiers to their associated
    sequencing groups (sample_to_sg_attribute_map).
    """
    logger.info(f'Transferring {len(samples)} samples')
    SequencingGroupAttributes = dict[Tuple[str, str, str], str]
    sample_to_sg_attribute_map: dict[Tuple[str, str], SequencingGroupAttributes] = {}
    old_sid_to_new_sid: dict[str, str] = {}
    for s in samples:
        exid_old_sid = (s['externalId'], s['id'])
        if exid_old_sid not in sample_to_sg_attribute_map:
            sample_to_sg_attribute_map[exid_old_sid] = {}
        for sg in s['sequencingGroups']:
            sample_to_sg_attribute_map[exid_old_sid][
                (sg['type'], sg['platform'], sg['technology'])
            ] = sg['id']

        sample_type = None if s['type'] == 'None' else s['type']
        existing_sid: str | None = None
        if existing_sample := get_existing_sample(existing_data, s['externalId']):
            existing_sid = existing_sample['id']

        existing_pid: int | None = None
        if s['participant']:
            existing_pid = upserted_participant_map[s['participant']['externalId']]

        sample_upsert = SampleUpsert(
            external_ids=s['externalIds'],
            type=sample_type or None,
            meta=(copy_files_in_dict(s['meta'], project) or {}),
            participant_id=existing_pid,
            sequencing_groups=upsert_sequencing_groups(s, existing_data, project),
            id=existing_sid,
        )

        logger.info(f'Processing sample {s["id"]}')
        logger.info('Creating test sample entry')
        new_sample_id = sapi.create_sample(
            project=target_project,
            sample_upsert=sample_upsert,
        )
        old_sid_to_new_sid[s['id']] = new_sample_id

    return old_sid_to_new_sid, sample_to_sg_attribute_map


def upsert_sequencing_groups(
    sample: dict, existing_data: dict, project: str
) -> list[SequencingGroupUpsert]:
    """Create SG Upsert Objects for a sample"""
    sgs_to_upsert: list[SequencingGroupUpsert] = []
    for sg in sample.get('sequencingGroups'):
        existing_sg = get_existing_sg(
            existing_data, sample.get('externalId'), sg.get('type')
        )
        existing_sgid = existing_sg.get('id') if existing_sg else None
        sg_upsert = SequencingGroupUpsert(
            id=existing_sgid,
            external_ids=sg.get('externalIds') or {},
            meta=sg.get('meta'),
            platform=sg.get('platform'),
            technology=sg.get('technology'),
            type=sg.get('type'),
            assays=upsert_assays(
                sg, existing_sgid, existing_data, sample.get('externalId'), project
            ),
        )
        sgs_to_upsert.append(sg_upsert)

    return sgs_to_upsert


def upsert_assays(
    sg: dict,
    existing_sgid: str | None,
    existing_data: dict,
    sample_external_id,
    project: str,
) -> list[AssayUpsert]:
    """Create Assay Upsert Objects for a sequencing group"""
    print(sg)
    assays_to_upsert: list[AssayUpsert] = []
    _existing_assay: dict[str, str] = {}
    for assay in sg.get('assays'):
        # Check if assay exists
        if existing_sgid:
            _existing_assay = get_existing_assay(
                existing_data, sample_external_id, existing_sgid, assay
            )
        existing_assay_id = _existing_assay.get('id') if _existing_assay else None
        assay_upsert = AssayUpsert(
            type=assay.get('type'),
            id=existing_assay_id,
            external_ids=assay.get('externalIds') or {},
            meta=copy_files_in_dict(assay.get('meta'), project),
        )

        assays_to_upsert.append(assay_upsert)
    return assays_to_upsert


def get_new_sg_id(
    old_sid: str,
    new_sg_attributes: tuple[str, str, str],
    old_sid_to_new_sid: dict[str, str],
    sample_to_sg_attribute_map: dict[tuple, dict[tuple, str]],
    new_sg_data: dict[str, dict[str, Any]],
):
    """
    Returns the new sequencing group id for a given sample id and sequencing group attributes.
    Args:
        old_sid (str): The sample id to search for.
        new_sg_attributes (tuple[str, str, str]): The attributes of the sequencing group to search for.
        old_sid_to_new_sid (dict[str, str]): A map from old sample ids to new sample ids.
        sample_to_sg_attribute_map (dict[tuple, dict[tuple, str]]): A map from (seid, sid) keys to a map of sequencing group attribute keys to old sequencing group ids.
        new_sg_data (dict[dict, Any]): The data containing the new samples and their sequencing groups.


        Example:
        old_sid = 'XPGAAAAA'
        old_sid_to_new_sid = {'XPGAAAAA': 'XPGBBBBB',
                              'XPGCCCCC': 'XPGDDDDD',
                              'XPGEEEEE': 'XPGFFFFF',
        sample_to_sg_attribute_map = {('HG00011', 'XPGAAAAA'): {('exome', 'Illumina', 'short-read'): 'CPGAAAAA',
                                      ('HG00011', 'XPGCCCCC'): {('exome', 'Illumina', 'short-read'): 'CPGBBBBB', # same participant, different sample to above
                                                                ('genome', 'Illumina', 'short-read'): 'CPGCCCCC'} # same participant, same sample, different sequencing group to above
                                      ('HG00022', 'XPGEEEEE'): {('exome', 'Illumina', 'short-read'): 'CPGDDDDD'} # different participant
        new_sg_data = {'project': {'samples': [{'id': 'XPGBBBBB', 'externalId': 'HG00011', 'sequencingGroups': [{'id': 'CPGZZZZZ', 'type': 'exome', 'platform': 'Illumina', 'technology': 'short-read'}],
                                               {'id': 'XPGDDDDD', 'externalId': 'HG00011', 'sequencingGroups': [{'id': 'CPGYYYYY', 'type': 'exome', 'platform': 'Illumina', 'technology': 'short-read'}]}}
                                               {'id': 'XPGDDDDD', 'externalId': 'HG00011', 'sequencingGroups': [{'id': 'CPGXXXXX', 'type': 'genome', 'platform': 'Illumina', 'technology': 'short-read'}]}}
                                               {'id': 'XPGFFFFF', 'externalId': 'HG00022', 'sequencingGroups': [{'id': 'CPGWWWWW', 'type': 'exome', 'platform': 'Illumina', 'technology': 'short-read'}]}}

        new_sg_attributes = ('exome', 'Illumina', 'short-read')

        The new_sg_attributes maps the newly created sequencing group to the ('externalId', 'old_sample_id') of the specific sample of the participant it was created from.

    """
    if old_sid in old_sid_to_new_sid:
        if not isinstance(new_sg_data.get('project', {}).get('samples'), list):
            raise TypeError("'samples' must be a list")
        for new_sample in new_sg_data['project']['samples']:
            if (
                new_sample['id'] == old_sid_to_new_sid[old_sid]
            ):  # If new sample maps to old sample
                for new_sg in new_sample['sequencingGroups']:
                    seid_sid_key = (new_sample['externalId'], old_sid)
                    sg_attribute_key = (
                        new_sg['type'],
                        new_sg['platform'],
                        new_sg['technology'],
                    )
                    if (
                        sample_to_sg_attribute_map[seid_sid_key].get(new_sg_attributes)
                        and new_sg_attributes == sg_attribute_key
                    ):
                        new_sg_id = new_sg['id']
                        return [new_sg_id]
                    # continue to next sg in new_sample
                    continue

    return ValueError(
        f'Old sample {old_sid} not found in mapping of old to new sample ids.'
    )


def transfer_analyses(
    samples: dict,
    existing_data: dict,
    target_project: str,
    project: str,
    old_sid_to_new_sid: dict[str, str],
    sample_to_sg_attribute_map: dict[tuple, dict[tuple, str]],
    update_embedded_ids: bool,
):
    """
    This function will transfer the analyses from the original project to the test project.
    """
    new_sg_data = query(SG_ID_QUERY, {'project': target_project})

    new_sg_map = {}
    for s in new_sg_data.get('project').get('samples'):
        sg_ids: list = []
        for sg in s.get('sequencingGroups'):
            sg_ids.append(sg.get('id'))
        new_sg_map[s.get('externalId')] = sg_ids

    for s in samples:
        for sg in s['sequencingGroups']:
            new_sg_attributes = sg.get('type'), sg.get('platform'), sg.get('technology')
            new_sequencing_group_id = get_new_sg_id(
                s['id'],
                new_sg_attributes,
                old_sid_to_new_sid,
                sample_to_sg_attribute_map,
                new_sg_data,
            )
            assert (
                len(new_sequencing_group_id) == 1
            ), f'Expected 1 new sequencing group id, got {len(new_sequencing_group_id)}'
            existing_sg = get_existing_sg(
                existing_data, s.get('externalId'), sg.get('type')
            )
            existing_sgid = existing_sg.get('id') if existing_sg else None
            for analysis in sg['analyses']:
                existing_analysis: dict = {}
                if existing_sgid:
                    existing_analysis = get_existing_analysis(
                        existing_data, s['externalId'], existing_sgid, analysis['type']
                    )
                existing_analysis_id = (
                    existing_analysis.get('id') if existing_analysis else None
                )
                if existing_analysis_id:
                    am = AnalysisUpdateModel(
                        type=analysis['type'],
                        output=copy_files_in_dict(
                            analysis['output'],
                            project,
                            (str(sg['id']), new_sg_map[s['externalId']][0]),
                            update_embedded_ids,
                        ),
                        status=AnalysisStatus(
                            analysis['status'].lower().replace('_', '-')
                        ),
                        sequencing_group_ids=new_sg_map[s['externalId']],
                        meta=analysis['meta'],
                    )
                    aapi.update_analysis(
                        analysis_id=existing_analysis_id,
                        analysis_update_model=am,
                    )
                else:
                    am = Analysis(
                        type=analysis['type'],
                        output=copy_files_in_dict(
                            analysis['output'],
                            project,
                            (str(sg['id']), new_sg_map[s['externalId']][0]),
                            update_embedded_ids,
                        ),
                        status=AnalysisStatus(
                            analysis['status'].lower().replace('_', '-')
                        ),
                        sequencing_group_ids=new_sequencing_group_id,
                        meta=analysis['meta'],
                    )

                    logger.info(f'Creating {analysis["type"]}analysis entry in test')
                    aapi.create_analysis(project=target_project, analysis=am)


def get_existing_sample(data: dict, sample_id: str) -> dict | None:
    """
    Get the existing sample object for this ID
    Returns:
        The Sample dictionary, or None if unmatched
    """
    for sample in data.get('project', {}).get('samples', []):
        if sample.get('externalId') == sample_id:
            return sample

    return None


def get_existing_sg(
    existing_data: dict, sample_id: str, sg_type: str = None, sg_id: str = None
) -> dict | None:
    """
    Find a SG ID in the main data based on a sample ID
    Match either on CPG ID or type (exome/genome)
    Returns:
        The SG Data, or None if no match is found
    """
    if not sg_type and not sg_id:
        raise ValueError('Must provide sg_type or sg_id when getting existing sg')
    if sample := get_existing_sample(existing_data, sample_id):
        for sg in sample.get('sequencingGroups'):
            if sg_id and sg.get('id') == sg_id:
                return sg
            if sg_type and sg.get('type') == sg_type:
                return sg

    return None


def get_existing_assay(
    data: dict, sample_id: str, sg_id: str, original_assay: dict
) -> dict | None:
    """
    Find assay in main data for this sample
    Returns:
        The Assay Data, or None if no match is found
    """
    if sg := get_existing_sg(existing_data=data, sample_id=sample_id, sg_id=sg_id):
        for assay in sg.get('assays', []):
            if assay.get('type') == original_assay.get('type') and assay.get(
                'meta'
            ) == original_assay.get('meta'):
                return assay

    return None


def get_existing_analysis(
    data: dict, sample_id: str, sg_id: str, analysis_type: str
) -> dict | None:
    """
    Find the existing SG for this sample, then identify any relevant analysis objs
    Returns:
        an analysis dict, or None if the right type isn't found
    """
    if sg := get_existing_sg(existing_data=data, sample_id=sample_id, sg_id=sg_id):
        for analysis in sg.get('analyses', []):
            if analysis.get('type') == analysis_type:
                return analysis
    return None


def get_sids_for_families(
    project: str, families_n: int, additional_families: set[str]
) -> set[str]:
    """
    Returns specific samples to be included in the test project.
    This process returns internal sample IDs for members of the families selected
    Remove any families with no samples from consideration prior to random selection
    """

    family_sgid_output = query(QUERY_FAMILY_SAMPLES, {'project': project})

    all_family_sgids = family_sgid_output.get('project', {}).get('families', [])
    assert all_family_sgids, 'No families returned in GQL result'

    # subselect down to families with samples
    families_with_samples = [
        fam for fam in all_family_sgids if any(participant.get('samples') for participant in fam['participants'])
    ]

    # 1. Remove the specifically requested families
    user_input_families = [
        fam for fam in families_with_samples if fam['externalId'] in additional_families
    ]

    # TODO: Replace this with the nice script that randomly selects better :)
    # 2. Randomly select from the remaining families (families_n can be 0)
    user_input_families.extend(
        random.sample(
            [
                fam
                for fam in families_with_samples
                if fam['externalId'] not in additional_families
            ],
            families_n,
        )
    )

    # log the gathered families, and any requested families we couldn't find
    selected_family_ids = [fam['externalId'] for fam in user_input_families]
    logger.info(f'Families selected: {", ".join(selected_family_ids)}')
    if failed_to_find := [user_fam_id for user_fam_id in additional_families if user_fam_id not in selected_family_ids]:
        logger.warning(f'Families not found: {", ".join(len(failed_to_find))}')

    # 3. Pull SGs from random + specific families
    included_sids: set[str] = set()
    for fam in user_input_families:
        for participant in fam['participants']:
            for sample in participant['samples']:
                included_sids.add(sample['id'])

    logger.info(f'Identified {len(included_sids)} samples from {len(user_input_families)} families')

    return included_sids


def get_sids_for_cohorts(
    project: str, cohorts: set[str], cohort_samples_n: int
) -> set[str]:
    """Returns cohort_samples_n specific samples for given cohort IDs."""

    cohort_sid_output = query(COHORT_QUERY, {'project': project})

    all_cohort_groups = cohort_sid_output.get('project', {}).get('cohorts', [])

    all_cohorts_sample_ids_subset: set[str] = set()
    for cohort in all_cohort_groups:
        sids_for_cohort: list[str] = []
        if cohort.get('id') in cohorts:
            seq_groups = cohort.get('sequencingGroups', [])
            for seq_group in seq_groups:
                sample = seq_group.get('sample')
                sids_for_cohort.append(sample['id'])
            all_cohorts_sample_ids_subset.update(
                random.sample(sids_for_cohort, cohort_samples_n)
            )

    return all_cohorts_sample_ids_subset


def transfer_families(
    initial_project: str, target_project: str, internal_participant_ids: list[int]
) -> list[int]:
    """Pull relevant families from the input project, and copy to target_project"""
    family_data = query(
        QUERY_FAMILY_BY_INTERNAL_IDS,
        variables={
            'project': initial_project,
            'internal_ids': internal_participant_ids,
        },
    )

    families: dict[int, dict] = {}
    for participant in family_data['project']['participants']:
        for family in participant['families']:
            families[family['id']] = family

    family_ids: list[int] = list(families.keys())

    tmp_family_tsv = 'tmp_families.tsv'
    # TODO: this doesn't match the default ordering in fapi.import_families, and is not passed as a list of headers
    family_tsv_headers = ['Family ID', 'Description', 'Coded Phenotype', 'Display Name']
    # Work-around as import_families takes a file.
    with open(tmp_family_tsv, 'wt') as tmp_families:
        tsv_writer = csv.writer(tmp_families, delimiter='\t')
        tsv_writer.writerow(family_tsv_headers)
        for family_id, family in families.items():
            # TODO only 3 elements are written to this 4-col TSV
            tsv_writer.writerow(
                [
                    # import_families() only imports the primary eid anyway
                    family['externalIds'][PRIMARY_EXTERNAL_ORG],
                    family['description'] or '',
                    family['codedPhenotype'] or '',
                ]
            )

    with open(tmp_family_tsv) as family_file:
        fapi.import_families(file=family_file, project=target_project)

    return family_ids


def transfer_ped(
    initial_project: str, target_project: str, family_ids: list[int]
) -> dict[str, int]:
    """Pull pedigree from the input project, and copy to target_project"""
    ped_tsv = fapi.get_pedigree(
        initial_project,
        export_type='tsv',
        internal_family_ids=family_ids,
    )
    tmp_ped_tsv = 'tmp_ped.tsv'
    # Work-around as import_pedigree takes a file.
    with open(tmp_ped_tsv, 'w') as tmp_ped:
        tmp_ped.write(ped_tsv)

    with open(tmp_ped_tsv) as ped_file:
        fapi.import_pedigree(
            file=ped_file,
            has_header=True,
            project=target_project,
            create_missing_participants=True,
        )

    # Get map of external participant id to internal
    participant_output = query(PARTICIPANT_QUERY, {'project': target_project})
    return {
        external_id: participant['id']
        for participant in participant_output.get('project').get('participants')
        for external_id in participant['externalIds'].values()
    }


def transfer_participants(
    target_project: str,
    participant_data,
) -> dict[str, int]:
    """Transfers relevant participants between projects"""

    existing_participants = query(PARTICIPANT_QUERY, {'project': target_project})

    target_project_pid_map = {
        external_id: participant['id']
        for participant in existing_participants.get('project').get('participants')
        for external_id in participant['externalIds'].values()
    }

    participants_to_transfer = []
    for participant in participant_data:
        for alt_id in participant['externalIds'].values():
            if alt_id in target_project_pid_map:
                # Participants with id field will be updated & those without will be inserted
                participant['id'] = target_project_pid_map[alt_id]
            else:
                del participant['id']
            transfer_participant = {
                'external_ids': participant['externalIds'],
                'meta': participant.get('meta') or {},
                'karyotype': participant.get('karyotype'),
                'reported_gender': participant.get('reportedGender'),
                'reported_sex': participant.get('reportedSex'),
                'id': participant.get('id'),
                'samples': [],
            }
            # Participants are being created before the samples are, so this will be empty for now.
            participants_to_transfer.append(transfer_participant)

    upserted_participants = papi.upsert_participants(
        target_project, participant_upsert=participants_to_transfer
    )

    external_to_internal_participant_id_map: dict[str, int] = {}

    for participant in upserted_participants:
        for external_id in participant['externalIds'].values():
            external_to_internal_participant_id_map[external_id] = participant['id']

    return external_to_internal_participant_id_map


def get_random_families(
    families: list[dict[str, str]],
    families_n: int,
    include_single_person_families: bool = False,
) -> list[str]:
    """Obtains a subset of families, that are a little less random.
    By default single-person families are discarded.
    The function aims to evenly distribute the families chosen by size.
    For example, if the composition of families inputted is as follows
    Duos - 5 families, Trios - 10 families, Quads - 5 families
    and families_n = 4
    Then this function will randomly select, 1 duo, 1 quad, and 2 trios.
    """

    family_sizes = dict(Counter([fam['family_id'] for fam in families]))

    # Discard single-person families
    family_threshold = 0 if include_single_person_families else 1
    families_within_threshold = [
        k for k, v in family_sizes.items() if v > family_threshold
    ]

    # Get family size distribution, i.e. {1:[FAM1, FAM2], 2:[FAM3], 3:[FAM4,FAM5, FAM6]}
    distributed_by_size: dict[int, list[str]] = {}
    for k, v in family_sizes.items():
        if k in families_within_threshold:
            if distributed_by_size.get(v):
                distributed_by_size[v].append(k)
            else:
                distributed_by_size[v] = [k]

    sizes = len(list(distributed_by_size.keys()))
    returned_families: list[str] = []

    proportion = families_n / len(families_within_threshold)

    if sizes <= families_n:
        for _s, fams in distributed_by_size.items():
            n_pull = round(proportion * len(fams))
            returned_families.extend(random.sample(fams, n_pull))

    else:
        # we can't evenly distribute, so we'll just pull randomly
        returned_families = random.sample(families_within_threshold, families_n)

    return returned_families


def file_reheaderable(path: str) -> bool:
    """Returns whether this is a file that can be reheadered"""
    return any(path.endswith(e) for e in ['.cram', '.vcf.gz', '.crai', '.tbi', '.md5'])


def main_file_commands(batch, job, old_path: str, new_path: str, sid: tuple[str, str]):
    """Adds the commands required to reheader the main CRAM/VCF file"""
    if new_path.endswith('.cram'):
        job.image(image_path('samtools'))
        old_file = batch.read_input(old_path)
        job.declare_resource_group(
            new_file={'main': '{root}.cram', 'crai': '{root}.cram.crai'}
        )
        job.command(rf"""
            samtools reheader --no-PG --in-place --command 'sed /^@RG/s/\<SM:{sid[0]}\>/SM:{sid[1]}/g' {old_file}
            mv {old_file} {job.new_file.main}
        """)
        batch.write_output(job.new_file.main, new_path)

    elif new_path.endswith('.vcf.gz'):
        job.image(image_path('bcftools'))
        old_file = batch.read_input(old_path)
        job.declare_resource_group(
            new_file={'main': '{root}.vcf.gz', 'tbi': '{root}.vcf.gz.tbi'}
        )
        job.command(rf"""
            echo {sid[0]} {sid[1]} > $BATCH_TMPDIR/samples.txt
            bcftools reheader --samples $BATCH_TMPDIR/samples.txt -o {job.new_file.main} {old_file}
        """)
        batch.write_output(job.new_file.main, new_path)


def extra_file_commands(batch, job, new_path: str):
    """Adds the commands required to regenerate the various associated files"""
    if new_path.endswith('.crai'):
        job.command(rf'samtools index -o {job.new_file.crai} {job.new_file.main}')
        batch.write_output(job.new_file.crai, new_path)

    elif new_path.endswith('.tbi'):
        job.command(rf"""
            tabix {job.new_file.main}
            # Mention {job.new_file.tbi} for hail as tabix has no -o option
        """)
        batch.write_output(job.new_file.tbi, new_path)

    elif new_path.endswith('.md5'):
        job.command(rf"md5sum {job.new_file.main} | cut -d' ' -f1 > {job.new_md5}")
        batch.write_output(job.new_md5, new_path)


def copy_files_in_dict(
    d,
    dataset: str,
    sid_replacement: tuple[str, str] = None,
    update_embedded_ids: bool = False,
):
    """
    Replaces all `gs://cpg-{project}-main*/` paths
    into `gs://cpg-{project}-test*/` and creates copies if needed
    If `d` is dict or list, recursively calls this function on every element
    If `d` is str, replaces the path
    """
    if not d:
        return d
    if isinstance(d, str) and d.startswith(f'gs://cpg-{dataset}-main'):
        logger.info(f'Looking for analysis file {d}')
        old_path = d
        if not file_exists(old_path):
            logger.warning(f'File {old_path} does not exist')
            return d
        new_path = old_path.replace(
            f'gs://cpg-{dataset}-main', f'gs://cpg-{dataset}-test'
        )
        # Replace the internal sample ID from the original project
        # With the new internal sample ID from the test project
        if sid_replacement is not None:
            new_path = new_path.replace(sid_replacement[0], sid_replacement[1])

        job = None

        if not file_exists(new_path):
            if update_embedded_ids and sid_replacement and file_reheaderable(new_path):
                logger.info(f'Adding job to rewrite {old_path} to {new_path}')
                batch = get_batch()
                job = batch.new_job(f'Reheader {new_path} from {old_path}')
                main_file_commands(batch, job, old_path, new_path, sid_replacement)
                job.storage(file_size(old_path) + 2 * 1024**3)  # Allow 2 GiB extra
            else:
                cmd = ['gcloud', 'storage', 'cp', old_path, new_path]
                logger.info(f'Copying file in metadata: {" ".join(cmd)}')
                subprocess.run(cmd, check=True)
        else:
            logger.error(f'{new_path} already exists')

        extra_exts = ['.md5']
        if new_path.endswith('.vcf.gz'):
            extra_exts.append('.tbi')
        if new_path.endswith('.cram'):
            extra_exts.append('.crai')
        for ext in extra_exts:
            if file_exists(old_path + ext) and not file_exists(new_path + ext):
                if job and file_reheaderable(new_path + ext):
                    logger.info(f'Extending job to regenerate {new_path + ext}')
                    extra_file_commands(batch, job, new_path + ext)
                else:
                    cmd = ['gcloud', 'storage', 'cp', old_path + ext, new_path + ext]
                    logger.info(f'Copying extra file in metadata: {" ".join(cmd)}')
                    subprocess.run(cmd, check=True)
            elif file_exists(old_path + ext):
                logger.error(f'{new_path} already exists')
        return new_path
    if isinstance(d, list):
        return [copy_files_in_dict(x, dataset) for x in d]
    if isinstance(d, dict):
        return {k: copy_files_in_dict(v, dataset) for k, v in d.items()}
    return d


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
    :param path: path to the file/directory/object
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


def file_size(path: str) -> int:
    """Return the size (in bytes) of the object, which is assumed to exist"""
    if path.startswith('gs://'):
        (bucket, path) = path.removeprefix('gs://').split('/', maxsplit=1)
        return storage.Client().get_bucket(bucket).get_blob(path).size

    return os.path.getsize(path)


if __name__ == '__main__':
    parser = ArgumentParser(description='Argument parser for subset generator')
    parser.add_argument(
        '--project', required=True, help='The sample-metadata project ($DATASET)'
    )
    parser.add_argument('-n', type=int, help='# Random Samples to copy', default=0)
    parser.add_argument('-f', type=int, help='# Random families to copy', default=0)
    parser.add_argument(
        '--nsamples-cohort',
        type=int,
        help='# Random samples to copy from each cohort',
        default=0,
    )
    # Flag to be used when there isn't available pedigree/family information.
    parser.add_argument(
        '--skip-ped',
        action='store_true',
        help='Skip transferring pedigree/family information',
    )
    parser.add_argument(
        '--families',
        nargs='+',
        help='Additional families to include.',
        type=str,
        default={},
    )
    parser.add_argument(
        '--samples',
        nargs='+',
        help='Additional samples to include.',
        type=str,
        default={},
    )
    parser.add_argument(
        '--cohorts',
        nargs='+',
        help='Cohorts to take random samples from.',
        type=str,
        default={},
    )
    parser.add_argument(
        '--noninteractive', action='store_true', help='Skip interactive confirmation'
    )
    parser.add_argument(
        '--update-embedded-ids',
        action=BooleanOptionalAction,
        help='Update IDs embedded within CRAM and VCF files (on by default)',
        default=True,
    )
    args, fail = parser.parse_known_args()
    if fail:
        parser.print_help()
        raise AttributeError(f'Invalid arguments: {fail}')

    main(
        project=args.project,
        samples_n=args.n,
        families_n=args.f,
        cohort_samples_n=args.nsamples_cohort,
        additional_samples=set(args.samples),
        additional_families=set(args.families),
        skip_ped=args.skip_ped,
        cohorts=set(args.cohorts),
        update_embedded_ids=args.update_embedded_ids,
    )
