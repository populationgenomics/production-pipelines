#!/usr/bin/env python3


"""
Aim is to take all previously batched samples and re-batch them to
get the best possible groupings of samples for GATK-SV

To do this we join every QC table from Evidence QC which has run so far,
filter to retain only active SG IDs, then re-batch cleanly.
"""


import json
import logging
from argparse import ArgumentParser

import pandas as pd

from cpg_utils import to_path
from cpg_workflows.jobs.sample_batching import batch_sgs
from metamist.graphql import gql, query

FIND_ACTIVE_SGS = gql(
    """
query FindActiveSGs($project: String!) {
  project(name: $project) {
    sequencingGroups(activeOnly: {eq: true}) {
      id
    }
  }
}
""",
)


def collect_all_sgids(projects: list[str]) -> list[str]:
    """
    Collect all active SG IDs for a project

    Args:
        projects (list[str]): all the projects to collect samples for

    Returns:
        a list of all unique active sample IDs
    """

    ongoing_collection: set[str] = set()

    for each_project in projects:
        response = query(FIND_ACTIVE_SGS, {'project': each_project})
        ongoing_collection.update([sg['id'] for sg in response['project']['sequencingGroups']])
    return sorted(ongoing_collection)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the QC tables', nargs='+', required=True)
    parser.add_argument('-p', help='Names of all relevant projects', nargs='+', required=True)
    parser.add_argument('-o', help='Where to write the output', required=True)
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(f'Unrecognised arguments: {unknown}')

    assert args.p and args.i
    logging.info(f'Collecting all active SG IDs for {args.p}')
    all_sg_ids = collect_all_sgids(args.p)

    dataframes = []
    for each in args.i:
        logging.info(f'Loading {each}')
        this_df = pd.read_csv(each, sep='\t', low_memory=False)
        this_df.columns = [x.replace('#', '') for x in this_df.columns]

        # filter to the active SGs we're interested in
        this_df = this_df.query('ID in @all_sg_ids')
        dataframes.append(this_df)

    one_big_df = pd.concat(dataframes).drop_duplicates()

    if len(one_big_df) == 0:
        raise ValueError('No samples found in the QC tables')

    # now make some batches
    batches = batch_sgs(one_big_df, min_batch_size=200, max_batch_size=300)

    with to_path(args.o).open('w') as f:
        f.write(json.dumps(batches, indent=4))
