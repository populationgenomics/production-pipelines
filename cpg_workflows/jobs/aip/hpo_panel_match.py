"""
version of prepare_aip_cohort designed to be run as part of a continuous pipeline

minimally reproduces
populationgenomics/automated-interpretation-pipeline/blob/main/helpers/prepare_aip_cohort.py
populationgenomics/automated-interpretation-pipeline/blob/main/helpers/hpo_panel_matching.py
- query for participants and HPOs
- query for panels and HPOs
- write a new file containing participant-panel matches
"""


import json
import logging
import re
import requests
from argparse import ArgumentParser
from collections import defaultdict

import networkx
from obonet import read_obo

from metamist.graphql import gql, query


HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:[0-9]+')
MAX_DEPTH: int = 3
PANELS_ENDPOINT = 'https://panelapp.agha.umccr.org/api/v1/panels/'


def get_json_response(url: str) -> dict:
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return

    Args:
        url (str): str URL to retrieve JSON format data from

    Returns:
        the JSON response from the endpoint
    """

    response = requests.get(url, headers={'Accept': 'application/json'}, timeout=60)
    response.raise_for_status()
    return response.json()


def get_panels(endpoint: str = PANELS_ENDPOINT) -> dict[str, set[int]]:
    """
    query panelapp, and collect panels by HPO term

    Args:
        endpoint (str): URL for panels

    Returns:
        dict: {HPO_Term: [panel_ID, panel_ID],}
    """

    hpo_dict = defaultdict(set)

    while endpoint:
        endpoint_data = get_json_response(endpoint)
        for panel in endpoint_data['results']:

            # can be split over multiple strings
            relevant_disorders = ' '.join(panel['relevant_disorders'] or [])
            for match in re.findall(HPO_RE, relevant_disorders):
                hpo_dict[match].add(int(panel['id']))

        # cycle through additional pages
        # why don't GEL make the panelapp API public...
        if endpoint_data['next']:
            endpoint = endpoint_data['next']
        else:
            break

    return dict(hpo_dict)


def get_participant_hpos(dataset: str) -> tuple[dict, set[str]]:
    """
    gql query, pull out family details and HPO terms
    may be a little overloaded at the moment
    Args:
        dataset (str): dataset name

    Returns:
        dict of per-participant details, and set of all HPO terms
    """

    query_string = gql(
        """
        query MyQuery($project: String!) {
            project(name: $project) {
                sequencingGroups {
                    sample {
                        participant {
                            phenotypes
                            externalId
                            families {
                                externalId
                            }
                        }
                    }
                    id
                }
            }
        }"""
    )

    result = query(query_string, variables={'project': dataset})
    hpo_dict = {}
    all_hpo: set[str] = set()
    for sg in result['project']['sequencingGroups']:
        hpos = set(
            HPO_RE.findall(sg['sample']['participant']['phenotypes'].get(HPO_KEY, ''))
        )
        all_hpo.update(hpos)
        hpo_dict[sg['id']] = dict(
            hpo_terms=hpos,
            family_id=sg['sample']['participant']['families'][0]['externalId'],
            external_id=sg['sample']['participant']['externalId'],
            panels={137},  # baseline panel is always mendeliome
        )
    return hpo_dict, all_hpo


def match_hpo_terms(
    panel_map: dict[str, set[int]],
    hpo_tree: networkx.MultiDiGraph,
    hpo_str: str,
    max_layer_delta: int = 3,
    layers_scanned: int = 0,
    selections: set | None = None,
) -> set[int]:
    """
    get a list of panels which are relevant for this HPO
    this includes a manual recursive implementation of the edge traversal
    main reason is to take superseded terms into account

    instead of just checking parent(s), we also check if a term is obsolete
    if so, we instead check each replacement term

    for live terms we recurse on all parents

    relevant usage guide:
    https://github.com/dhimmel/obonet/blob/main/examples/go-obonet.ipynb
    """

    if selections is None:
        selections = set()

    # identify identical match and select the panel
    if hpo_str in panel_map:
        selections.update(panel_map[hpo_str])

    if layers_scanned >= max_layer_delta:
        return selections

    # if a node is invalid, recursively call this method for each replacement D:
    # there are simpler ways, just none that are as fun to write
    if not hpo_tree.has_node(hpo_str):
        logging.error(f'HPO term was absent from the tree: {hpo_str}')
        return selections

    hpo_node = hpo_tree.nodes[hpo_str]
    if hpo_node.get('is_obsolete', 'false') == 'true':
        for hpo_term in hpo_node.get('replaced_by', []):
            selections.update(
                match_hpo_terms(
                    panel_map,
                    hpo_tree,
                    hpo_term,
                    max_layer_delta,
                    layers_scanned + 1,
                    selections,
                )
            )
    # search for parent(s), even if the term is obsolete
    for hpo_term in hpo_node.get('is_a', []):
        selections.update(
            match_hpo_terms(
                panel_map,
                hpo_tree,
                hpo_term,
                max_layer_delta,
                layers_scanned + 1,
                selections,
            )
        )
    return selections


def match_hpos_to_panels(
    hpo_to_panel_map: dict, hpo_file: str, all_hpos: set[str]
) -> dict[str, set[int]]:
    """
    take the HPO terms from the participant metadata, and match to panels
    Args:
        hpo_to_panel_map (dict): panel IDs to all related panels
        hpo_file (str): path to an obo file containing HPO tree
        all_hpos (set[str]): collection of all unique HPO terms

    Returns:
        a dictionary linking all HPO terms to a corresponding set of Panel IDs
    """
    hpo_graph = read_obo(hpo_file, ignore_obsolete=False)

    hpo_to_panels = {}
    for hpo in all_hpos:
        panel_ids = match_hpo_terms(
            panel_map=hpo_to_panel_map,
            hpo_tree=hpo_graph,
            hpo_str=hpo,
            max_layer_delta=MAX_DEPTH,
        )
        hpo_to_panels[hpo] = panel_ids

    return hpo_to_panels


def match_participants_to_panels(participant_hpos: dict, hpo_panels: dict):
    """
    take the two maps of Participants: HPOs, and HPO: Panels
    blend the two to find panels per participant

    For each participant, find any HPO terms which were matched to panels
    for each matched term, add the panel(s) to the participant's private set

    Args:
        participant_hpos (dict): CPG ID to other details
        hpo_panels (dict): lookup of panels per HPO term
    """

    for party_data in participant_hpos.values():
        for hpo_term in party_data['hpo_terms']:
            if hpo_term in hpo_panels:
                party_data['panels'].update(hpo_panels[hpo_term])
        party_data['panels'] = list(party_data['panels'])
        party_data['hpo_terms'] = list(party_data['hpo_terms'])


def main(dataset: str, hpo_file: str, panel_out: str):
    """
    main method to do the fun things

    Args:
        dataset (str): name of the dataset to query
        hpo_file (str): path to a localised HPO OBO file
        panel_out (str): where to write final panel file
    """
    panels_by_hpo = get_panels()
    hpo_dict, all_hpo = get_participant_hpos(dataset=dataset)
    hpo_to_panels = match_hpos_to_panels(
        hpo_to_panel_map=panels_by_hpo, hpo_file=hpo_file, all_hpos=all_hpo
    )
    match_participants_to_panels(hpo_dict, hpo_to_panels)
    with open(panel_out, 'w', encoding='utf-8') as handle:
        json.dump(hpo_dict, handle, indent=4)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dataset', help='metamist project name')
    parser.add_argument('--hpo', help='local copy of HPO obo file')
    parser.add_argument('--out', help='panel file to write')
    args = parser.parse_args()
    main(dataset=args.dataset, hpo_file=args.hpo, panel_out=args.out)
