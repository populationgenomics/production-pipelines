#!/usr/bin/env python3

"""
Register reheadered NAGIM CRAM files into Metamist

Typical usage:

analysis-runner ... scripts/register_reheadered_nagim_crams.py
  --gcs-prefix gs://.../reheadered/
  --dry-run (optional)
"""

import datetime
import logging
import subprocess
from typing import Any

import click
import coloredlogs

from cpg_utils import to_path
from cpg_utils.config import get_config
from metamist import models
from metamist.apis import AnalysisApi
from metamist.graphql import gql, query

FIND_ANALYSES = gql(
    """
query AnalysesQuery($project: String!) {
  project(name: $project) {
    analyses(type: {eq: "cram"}, active: {eq: true}) {
      id
      meta
      output
      sequencingGroups {
        sample {
          id
        }
        id
      }
      status
      type
    }
  }
}
""",
)

FIND_SAMPLES = gql(
    """
query SamplesQuery($project: String!) {
  project(name: $project) {
    samples {
      id
      sequencingGroups(activeOnly: {eq: true}) {
        id
        analyses(type: {eq: "cram"}, active: {eq: true}) {
          id
          output
          status
        }
      }
    }
  }
}
""",
)

MUTATION_QUERY = gql(
    """
    mutation MyMutation($Content: String!, $Id: String!) {
      sequencingGroup {
        addComment(content: $Content, id: $Id) {
          content
          id
        }
      }
    }
    """,
)


def do_metamist_update(
    dry_run: bool,
    dataset: str,
    activeseqgroup: str,
    oldanalysis: dict[str, Any],
    oldpath: str,
    newpath: str,
):
    aapi = AnalysisApi()

    newmeta: dict[str, str] = oldanalysis['meta']
    newmeta['source'] = f'Reheadered NAGIM from {oldpath}'

    newanalysis = models.Analysis(
        type=oldanalysis['type'],
        status=models.AnalysisStatus(oldanalysis['status'].lower()),
        output=str(newpath),
        sequencing_group_ids=[activeseqgroup],
        meta=newmeta,
    )
    if dry_run:
        logging.info('Going to register the following new analysis:')
        logging.info(f'{newanalysis}')
        logging.info(f'Will add the following comment to the sequencing group: {activeseqgroup}:')
        logging.info(
            f'The reheadered NAGIM CRAM {newpath} has been registered as the latest analysis on {datetime.datetime.today()}, replacing the older Dragmap CRAM. A new Dragen 3_7_8 CRAM will replace this at a later date.',
        )
    else:
        aid = aapi.create_analysis(dataset, newanalysis)
        logging.info(f'Created Analysis(id={aid}, output={newpath}) in {dataset}')

        reheadered_comment = f'The reheadered NAGIM CRAM {newpath} has been registered as the latest analysis on {datetime.datetime.today()}, replacing the older Dragmap CRAM. A new Dragen 3_7_8 CRAM will replace this at a later date.'
        mutation_result = query(MUTATION_QUERY, variables={'Content': reheadered_comment, 'Id': activeseqgroup})


@click.command()
@click.option('--gcs-prefix', help='Prefix in GCP to the files to be registered.')
@click.option('--dry-run', is_flag=True, help='Display information only without making changes')
def main(
    gcs_prefix: str,
    dry_run: bool,
):
    # config = get_config(True)

    # dataset = config['workflow']['dataset']
    # if config['workflow']['access_level'] == 'test' and not dataset.endswith('-test'):
    #     dataset = f'{dataset}-test'
    coloredlogs.install(level=logging.INFO)

    dataset = 'tob-wgs-test'

    project_analyses = query(FIND_ANALYSES, {'project': dataset})['project']['analyses']
    analysis_by_path = {analysis['output']: analysis for analysis in project_analyses if 'nagim' in analysis['output']}
    project_samples = query(FIND_SAMPLES, {'project': dataset})['project']['samples']
    sample_by_id = {sample['id']: sample['sequencingGroups'][0]['id'] for sample in project_samples}

    cram_list = analysis_by_path.keys()

    flist = (
        subprocess.run(['gcloud', 'storage', 'ls', f'{gcs_prefix}*.cram'], capture_output=True)
        .stdout.decode(
            'utf-8',
        )
        .strip()
        .split('\n')
    )

    logging.info(f'Total reheadered NAGIM CRAM file entries found: {len(cram_list)}')

    for fname in cram_list:
        path = to_path(fname)
        reheadered_path = path.parent / 'reheadered' / path.name
        if str(reheadered_path) not in flist:
            continue
        try:
            real_sequnging_group: str = sample_by_id[analysis_by_path[fname]['sequencingGroups'][0]['sample']['id']]
        except KeyError:
            logging.info(f'The sample {analysis_by_path[fname]["sequencingGroups"][0]["sample"]["id"]} is not active.')
            continue
        do_metamist_update(
            dry_run=dry_run,
            dataset=dataset,
            activeseqgroup=real_sequnging_group,
            oldanalysis=analysis_by_path[fname],
            oldpath=fname,
            newpath=str(reheadered_path),
        )


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter  # click will add the arguments
