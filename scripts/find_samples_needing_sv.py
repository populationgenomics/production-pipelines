#!/usr/bin/python3

"""
a script to find samples that need to be run through the structural variant pipeline
(i.e. are not present yet in any extant callsets)

usage:
    python3 find_samples_needing_sv.py
    -n <max_#_samples>
    -o <output_path>
    <projects in order of preference>

The integer argument -n will be the maximum number of samples to return
If all eligible SG IDs in all listed projects are not enough to reach this number,
the samples that were successfully collected will be returned

optional arguments:
    -e <path to json list containing sequencing groups to exclude from search>
    -a <path to json list containing sequencing groups to include specifically>
"""
import re
from argparse import ArgumentParser
from random import sample

import json

from metamist.graphql import query, gql


FIND_ANALYSES = gql(
    """
query MyQuery ($project: String!, $type: String!) {
  project(name: $project) {
    analyses(
      status: {eq: COMPLETED}
      type: {eq: $type}
      active: {eq: true}
    ) {
      meta
      sequencingGroups {
        id
      }
    }
  }
}
"""
)

GENOME_SGS = gql(
    """
query SG_Query($project: String!) {
  project(name: $project) {
    sequencingGroups(activeOnly: {eq: true}, type: {eq: "genome"}) {
      id
      type
    }
  }
}"""
)

# Sequencing groups matching these patterns have consistently failed
# so they are being screened permanently
BORKED_SGS = ['CPG1', 'CPG5']


def get_all_sgs(project: str, exclude_sgs: set[str] | None = None) -> set[str]:
    """
    find all sgs in the project

    Args:
        project ():
        exclude_sgs ():

    Returns:

    """

    if exclude_sgs is None:
        exclude_sgs = set()

    sgs = set()

    for sg in query(GENOME_SGS, {'project': project})['project']['sequencingGroups']:
        if sg['id'] in exclude_sgs:
            continue
        if any(sg['id'].startswith(bork) for bork in BORKED_SGS):
            continue
        sgs.add(sg['id'])
    return sgs


def get_called_sgs(project: str) -> set[str]:
    """
    find all sequencing groups in the project with SV calls

    Args:
        project (str): the project ID to use

    Returns:
        a set of all sequencing groups in the project
    """

    sgs = set()

    for analysis in query(FIND_ANALYSES, {'project': project, 'type': 'sv'})['project']['analyses']:
        if analysis['meta'].get('stage') != 'AnnotateDatasetSv':
            continue

        sgs.update(analysis['meta']['sequencing_groups'])

    return sgs


def get_project_crams(project: str) -> set[str]:
    """
    find all SGs with crams in the project

    Args:
        project (str): the project ID to use

    Returns:
        a set of all crams in the project
    """

    crams = set()

    for analysis in query(FIND_ANALYSES, {'project': project, 'type': 'cram'})['project']['analyses']:
        if analysis['sequencingGroups']:
            crams.add(analysis['sequencingGroups'][0]['id'])
    return crams


def main(
    max_samples: int,
    output_path: str,
    projects: list[str],
    exclude: str = None,
    additional: str | None = None,
):
    """

    Args:
        max_samples ():
        output_path ():
        projects ():
        exclude (str | None): a json file containing sequencing groups to exclude from search
        additional (str | None): a json file containing sequencing groups to include specifically
    """

    collected_sgs: set[str] = set()

    # get excluded SGs if that's valid
    if exclude is not None:
        with open(exclude) as f:
            exclude_sgs = set(json.load(f))
    else:
        exclude_sgs = set()

    # iterate over projects in order
    for project in projects:

        # decide if its time to stop
        if len(collected_sgs) >= max_samples:
            break

        # find all SGs with a registered ready CRAM
        project_crams = get_project_crams(project)

        # find all SGs, with optional manual exclusions
        all_eligible_sgs = get_all_sgs(project=project, exclude_sgs=exclude_sgs)

        # find all SGs with SV calls
        called_sgs = get_called_sgs(project=project)

        # find all SGs with CRAMs & without SV calls
        project_samples_to_call = project_crams.intersection(all_eligible_sgs) - called_sgs
        print(f'{project}: {len(project_samples_to_call)} samples to call')

        if len(project_samples_to_call) + len(collected_sgs) > max_samples:
            more_samples = max_samples - len(collected_sgs)
            collected_sgs.update(sample(list(project_samples_to_call), more_samples))
            print(f'Only added {more_samples} samples from {project}')

        else:
            collected_sgs.update(project_samples_to_call)

    if additional is not None:
        with open(additional) as f:
            collected_sgs.update(json.load(f))

    # write out the results
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(list(collected_sgs), f, indent=4)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--min_samples', type=int, default=200)
    parser.add_argument('-o', '--output_path', default='samples_needing_sv.json')
    parser.add_argument('-e', '--exclude', required=False, default=None)
    parser.add_argument('-a', '--additional', required=False, default=None)
    parser.add_argument('projects', nargs='+')
    args = parser.parse_args()
    main(
        args.min_samples,
        args.output_path,
        args.projects,
        exclude=args.exclude,
        additional=args.additional,
    )
