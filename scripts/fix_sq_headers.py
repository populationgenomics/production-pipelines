#!/usr/bin/env python3
"""
Usage: fix_sq_headers list|update

Uses the usual config entries (workflow.dataset/.access_level)
"""

import hashlib
import os
import subprocess as sb
import sys

from cpg_utils.hail_batch import get_config
from metamist.graphql import query, gql

FIND_CRAMS = gql("""
query CramQuery($project: String!) {
  project(name: $project) {
    sequencingGroups {
      id
      analyses(type: {eq: "cram"}) {
        id
        meta
        output
        active
        type
      }
    }
  }
}
""")


KNOWN_REFS = {
    '49f7a7583b832aed0104b1901f9cf5fd': 'Homo_sapiens_assembly38_masked',
    '760a471d32c9a356ffb37f632bebfea3': 'Homo_sapiens_assembly38[unmasked]',
    }


def set_htslib_token():
    gcloud_cmd = ['gcloud', 'auth', 'application-default', 'print-access-token']
    token = sb.run(gcloud_cmd, check=True, capture_output=True, text=True).stdout
    os.environ['GCS_OAUTH_TOKEN'] = token


def parse_one_SQ(line: str):
    table = {}
    for field in line.strip().split('\t')[1:]:
        key, value = field.split(':', maxsplit=1)
        table[key] = value
    return table


def classify_SQ(path: str):
    """Reads the file's headers and recognises well-known reference collections"""
    run = sb.run(['samtools', 'head', path], capture_output=True, text=True)
    if (run.returncode != 0):
        return None  # Probably file not found

    m5 = {}
    for line in run.stdout.splitlines():
        if line.startswith('@SQ'):
            field = parse_one_SQ(line)
            m5[field['SN']] = field['M5']

    text = '/'.join([f'{name}:{m5[name]}' for name in sorted(m5.keys())])
    hashed = hashlib.md5(text.encode('ASCII')).hexdigest()

    return KNOWN_REFS.get(hashed, hashed)


if __name__ == '__main__':
    config = get_config(True)

    dataset = config['workflow']['dataset']
    if config['workflow']['access_level'] == 'test' and not dataset.endswith('-test'):
        dataset = f'{dataset}-test'

    print(f'Finding CRAMs for {dataset}')

    cram_list = []

    for seqgroup in query(FIND_CRAMS, {'project': dataset})['project']['sequencingGroups']:
        for analysis in seqgroup['analyses']:
            fname = analysis['output']
            if fname.endswith('.cram'):
                cram_list.append(fname)
            else:
                sgid = seqgroup['id']
                anid = analysis['id']
                print(f'Sequencing group {sgid} analysis {anid} is not CRAM: {fname}')

    print(f'Total CRAM file database entries found: {len(cram_list)}')

    command = sys.argv[1] if len(sys.argv) > 1 else 'help'
    if command == 'list':
        set_htslib_token()
        for path in cram_list:
            ref = classify_SQ(path) or '(file missing)'
            print(f'{path}: {ref}')

    else:
        sys.exit(__doc__)
