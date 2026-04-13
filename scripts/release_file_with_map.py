#!/usr/bin/env python3

"""
Script to take a VCF and a project
- gets the project ID mapping (CPG ID to external participant ID)
- writes the mapping to a release bucket
- copies the VCF and corresponding index to the same dir in the same release bucket
"""

import json
import os
import subprocess
from argparse import ArgumentParser
from datetime import datetime, timezone

from google.cloud import storage

from cpg_utils import to_path
from metamist.graphql import gql, query

TODAY = datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')

ID_QUERY = gql(
    """
query MyQuery($dataset: String!) {
    project(name: $dataset) {
    sequencingGroups {
      id
      sample {
        participant {
          externalId
        }
      }
    }
  }
}""",
)


def upload_metadata_to_release(dataset: str, billing_project: str | None):
    """Uploads the compressed metadata.zip into a directory with today's date in the release bucket"""

    results = query(ID_QUERY, {'dataset': dataset})['project']['sequencingGroups']
    id_map = {sg['id']: sg['sample']['participant']['externalId'] for sg in results}

    with open('id_mapping.json', 'w') as f:
        json.dump(id_map, f, indent=4)

    upload_path = os.path.join(TODAY, 'id_mapping.json')

    release_bucket = f'cpg-{dataset}-release'

    if billing_project:
        client = storage.Client(project=billing_project)
    else:
        client = storage.Client()

    bucket = client.bucket(release_bucket, user_project=billing_project)

    zip_blob = bucket.blob(upload_path)

    zip_blob.upload_from_filename('id_mapping.json')

    print(f'Uploaded SG mapping to gs://{os.path.join(release_bucket, upload_path)}')


def copy_vcf_to_release(vcf: str, dataset: str, billing_project: str):
    """Copies a specified vcf to the release bucket."""

    release_path = f'gs://cpg-{dataset}-release/{TODAY}/'

    vcf_name = to_path(vcf)

    for ext in ['', '.tbi']:

        subprocess.run(
            [  # noqa: S603, S607
                'gcloud',
                'storage',
                '--billing-project',
                billing_project,
                'cp',
                f'{vcf}{ext}',
                os.path.join(release_path, f'{vcf_name}{ext}'),
            ],
            check=True,
        )


def main(
    dataset: str,
    billing_project: str | None,
    vcf: str,
):
    """Creates the metadata files and saves them to the output path"""

    upload_metadata_to_release(dataset=dataset, billing_project=billing_project or dataset)
    copy_vcf_to_release(vcf, dataset, billing_project or dataset)


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--vcf', required=True)
    args = parser.parse_args()

    main(dataset=args.dataset, billing_project=args.billing_project, vcf=args.vcf)
