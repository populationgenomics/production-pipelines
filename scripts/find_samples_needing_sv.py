#!/usr/bin/python3

"""
a script to find samples that need to be run through the structural variant pipeline
(i.e. are not present yet in any extant callsets)

usage:
    python3 find_samples_needing_sv.py
    -n <min_#_samples>
    -o <output_path>
    <projects in order of preference>

optional argument:
    -e <path to json file containing sequencing groups to exclude from search>
    -a <path to json file containing sequencing groups to include specifically>
"""
from argparse import ArgumentParser
from random import sample

import json

from metamist.graphql import query, gql


FIND_ANALYSES = gql(
    """
query MyQuery ($project: String!) {
  project(name: $project) {
    analyses(
      status: {eq: COMPLETED}
      type: {eq: "sv"}
      active: {eq: true}
    ) {
      meta
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
        if sg['id'] not in exclude_sgs:
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

    for analysis in query(FIND_ANALYSES, {'project': project})['project']['analyses']:
        if analysis['meta'].get('type') != 'annotated-sv-dataset-callset':
            continue
        assert analysis['meta'].get('stage') == 'AnnotateDatasetSv'

        sgs.update(analysis['meta']['sequencing_groups'])
    return sgs


def main(
    min_samples: int,
    output_path: str,
    projects: list[str],
    exclude: str = None,
    additional: str | None = None,
):
    """

    Args:
        min_samples ():
        output_path ():
        projects ():
        exclude (str | None): a json file containing sequencing groups to exclude from search
        additional (str | None): a json file containing sequencing groups to include specifically
    """

    collected_sgs = set()

    # get excluded SGs if that's valid
    if exclude is not None:
        with open(exclude) as f:
            exclude_sgs = set(json.load(f))
    else:
        exclude_sgs = set()

    # iterate over projects in order
    for project in projects:

        if len(collected_sgs) >= min_samples:
            break

        # find all SGs
        all_eligible_sgs = get_all_sgs(project=project, exclude_sgs=exclude_sgs)

        # find all SGs with SV calls
        called_sgs = get_called_sgs(project=project)

        # find all SGs without SV calls
        project_samples_to_call = all_eligible_sgs - called_sgs
        print(f'{project}: {len(project_samples_to_call)} samples to call')

        if len(project_samples_to_call) + len(collected_sgs) > min_samples:
            more_samples = min_samples - len(collected_sgs)
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
