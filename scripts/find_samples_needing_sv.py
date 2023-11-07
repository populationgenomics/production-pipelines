#!/usr/bin/python3

"""
a script to find samples that need to be run through the structural variant pipeline
(i.e. are not present yet in any extant callsets)

usage:
    python3 find_samples_needing_sv.py
    -n <min_#_samples>
    -o <output_path>
    <projects in order of preference>
"""

from argparse import ArgumentParser
import json
from metamist.graphql import query, gql


QUERY_STRING = f"""
"""


def main(min_samples: int, output_path: str, projects: list[str]):
    """

    Args:
        min_samples ():
        output_path ():
        projects ():

    Returns:

    """


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--min_samples', type=int, default=200)
    parser.add_argument('-o', '--output_path', default='samples_needing_sv.json')
    parser.add_argument('projects', nargs='+')
    args = parser.parse_args()
    main(args.min_samples, args.output_path, args.projects)
