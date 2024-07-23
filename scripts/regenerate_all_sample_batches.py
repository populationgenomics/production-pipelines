#!/usr/bin/env python3


"""
Aim is to take all previously batched samples and re-batch them to
get the best possible groupings of samples for GATK-SV

To do this we join every QC table from Evidence QC which has run so far,
filter to retain only active SG IDs, then re-batch cleanly.
"""


import json
from argparse import ArgumentParser

import pandas as pd

from cpg_utils import to_path
from cpg_workflows.jobs.sample_batching import batch_sgs, batch_sgs_by_library
from cpg_workflows.utils import get_logger
from metamist.apis import FamilyApi
from metamist.graphql import gql, query

FIND_ACTIVE_SGS = gql(
    """
query FindActiveSGs($project: String!) {
  project(name: $project) {
    sequencingGroups(activeOnly: {eq: true}) {
      id
      meta
    }
  }
}
""",
)

SG_PARTICIPANTS_QUERY = gql(
    """
    query FindSGParticipants($project: String!) {
        project(name: $project) {
            sequencingGroups {
                id
                sample {
                    participant {
                        id
                        externalId
                    }
                }
            }
        }
    }
    """,
)


def parse_sg_meta(df: pd.DataFrame, sg_meta: dict[str, dict[str, str]]) -> pd.DataFrame:
    """
    Parse the SG meta fields (sequencing facility, library) into the QC table

    Args:
        df (pd.DataFrame): the QC table
        sg_meta (dict[str, dict[str, str]]): the metadata for each SG

    Returns:
        the QC table with the meta fields added
    """
    # library might be "library_type" or "sequencing_library"
    # facility might be "facility" or "sequencing_facility"
    df['library'] = df['ID'].map(
        lambda x: sg_meta[x].get('library_type', sg_meta[x].get('sequencing_library', 'unknown')),
    )
    df['facility'] = df['ID'].map(
        lambda x: sg_meta[x].get('facility', sg_meta[x].get('sequencing_facility', 'unknown')),
    )
    return df


def get_sg_participants(project: str) -> list[dict]:
    """
    Get the participant ID for each SG using GraphQL

    Args:
        project (str): the project to collect samples for

    Returns:
        a list of dictionaries with the SG ID and the participant ID
    """
    response = query(SG_PARTICIPANTS_QUERY, {'project': project})
    return [
        {
            sg['id']: sg['sample']['participant']['externalId'],
        }
        for sg in response['project']['sequencingGroups']
    ]


def get_sg_families(projects: list[str]) -> dict[str, int]:
    """
    Get the family ID for each SG, by first getting the participant ID, then getting the
    pedigree for the project and mapping the participant ID to the family ID

    Args:
        projects (list[str]): all the projects

    Returns:
        a dictionary of SG ID to internal family ID
    """
    sg_families: dict[str, int] = {}
    participant_families: dict[str, int] = {}
    sg_participants: dict[str, str] = {}
    family_api = FamilyApi()
    for project in projects:
        sg_participants.update({x['id']: x['externalId'] for x in get_sg_participants(project)})
        participant_families.update(
            {
                x['individual_id']: int(x['family_id'])
                for x in family_api.get_pedigree(project=project, replace_with_family_external_ids=False)
            },
        )

    for sg, participant in sg_participants.items():
        sg_families[sg] = participant_families.get(participant, -1)
    return sg_families


def collect_all_sgids(projects: list[str]) -> tuple[list[str], dict[str, dict[str, str]]]:
    """
    Collect all active SG IDs for a project

    Args:
        projects (list[str]): all the projects to collect samples for

    Returns:
        a list of all unique active sample IDs
    """
    sg_meta: dict[str, dict[str, str]] = {}
    ongoing_collection: set[str] = set()

    for each_project in projects:
        response = query(FIND_ACTIVE_SGS, {'project': each_project})
        ongoing_collection.update([sg['id'] for sg in response['project']['sequencingGroups']])
        sg_meta.update({sg['id']: sg['meta'] for sg in response['project']['sequencingGroups']})
    return sorted(ongoing_collection), sg_meta


if __name__ == '__main__':
    get_logger(__file__).info('Starting the re-batching process')
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the QC tables', nargs='+', required=True)
    parser.add_argument('-p', help='Names of all relevant projects', nargs='+', required=True)
    parser.add_argument('-o', help='Where to write the output', required=True)
    parser.add_argument('--meta', help='Add the sg meta fields to the output', default=False, action='store_true')
    parser.add_argument('--family', help='Get family ID for each SG', default=False, action='store_true')
    parser.add_argument('--min', help='Min Batch Size', type=int, default=200)
    parser.add_argument('--max', help='Max Batch Size', type=int, default=300)
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(f'Unrecognised arguments: {unknown}')

    assert args.p and args.i
    get_logger().info(f'Collecting all active SG IDs for {args.p}')
    all_sg_ids, all_sg_meta = collect_all_sgids(args.p)
    get_logger().info(f'Identified {len(all_sg_ids)} active SG IDs')
    if args.family:
        get_logger().info('Collecting family IDs for each SG')
        sg_families = get_sg_families(args.p)

    dataframes = []
    for each in args.i:
        get_logger().info(f'Loading {each}')
        this_df = pd.read_csv(each, sep='\t', low_memory=False)
        this_df.columns = [x.replace('#', '') for x in this_df.columns]  # type: ignore

        # filter to the active SGs we're interested in
        this_df = this_df.query('ID in @all_sg_ids')
        if args.meta:
            this_df = parse_sg_meta(this_df, all_sg_meta)
        if args.family:
            this_df['family'] = this_df['ID'].map(lambda x: sg_families[x])
        dataframes.append(this_df)

    one_big_df = pd.concat(dataframes).drop_duplicates()

    if len(one_big_df) == 0:
        raise ValueError('No samples found in the QC tables')

    # now make some batches
    if args.meta:
        batches = batch_sgs_by_library(one_big_df, min_batch_size=args.min, max_batch_size=args.max)
    else:
        batches = batch_sgs(one_big_df, min_batch_size=args.min, max_batch_size=args.max)

    with to_path(args.o).open('w') as f:
        f.write(json.dumps(batches, indent=4))
