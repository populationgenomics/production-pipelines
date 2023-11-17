#!/usr/bin/env python3
"""
Usage: fix_sq_headers list|update

Uses the usual config entries (workflow.dataset/.access_level) and:

  images.samtools      Path of an image containing samtools
  references.expected  Name (from KNOWN_REFS) of the (bad) refset expected
  references.new       Dict file containing the replacement @SQ headers
"""

import hashlib
import os
import subprocess as sb
import sys

from hailtop.batch.resource import JobResourceFile

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import get_batch, get_config, image_path
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


def path_exists(path: Path):
    try:
        return path.exists()
    except Exception:
        # GSPath.exists() raises google.api_core.exceptions.NotFound if the bucket is missing
        return False


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


def do_reheader(newref, incram, outcram, outcrai):
    translate = {}
    with open(newref) as refp:
        for line in refp:
            if line.startswith('@SQ'):
                field = parse_one_SQ(line)
                translate[field['SN']] = (line.strip(), field['LN'])

    head_cmd = ['samtools', 'head', incram]
    headers = sb.run(head_cmd, check=True, capture_output=True, text=True).stdout

    newhdr_file = 'tmp.hdr'
    with open(newhdr_file, 'w') as newhdr:
        for line in headers.splitlines():
            if line.startswith('@SQ'):
                field = parse_one_SQ(line)
                new_line, new_len = translate[field['SN']]
                if field['LN'] != new_len:
                    raise ValueError(f'{incram}:{field["SN"]}: new length {new_len} differs '
                                     f'from existing {field["LN"]}')
                print(new_line, file=newhdr)
            else:
                print(line, file=newhdr)

    sb.run(['samtools', 'reheader', '--no-PG', '--in-place', newhdr_file, incram], check=True)
    os.rename(incram, outcram)
    sb.run(['samtools', 'index', '-o', outcrai, outcram], check=True)


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

    elif command == 'update':
        expected_refset = config['references']['expected']
        new_reference_path = config['references']['new']

        set_htslib_token()
        reheader_list = []
        for fname in cram_list:
            path = to_path(fname)
            reheadered_path = path.parent / 'reheadered' / path.name

            if not path_exists(path):
                print(f'{fname}: file does not exist')
            elif path_exists(reheadered_path):
                print(f'{reheadered_path}: reheadered file already exists')
            elif (ref := classify_SQ(fname)) != expected_refset:
                print(f'{fname}: unexpected original reference headers: {ref}')
            else:
                reheader_list.append((path, reheadered_path, path.stat().st_size))

        if len(reheader_list) == 0:
            sys.exit('No CRAMs to be reheadered!')

        print(f'CRAM files to be reheadered: {len(reheader_list)}')

        b = get_batch(f'Reheader {dataset} to {new_reference_path}')
        new_reference = b.read_input(new_reference_path)

        for path, newpath, filesize in reheader_list:
            j = b.new_python_job(f'Reheader {path}')
            j.image(image_path('samtools'))
            j.storage(filesize * 1.2)

            j.call(do_reheader, new_reference, b.read_input(str(path)), j.out_cram, j.out_crai)

            assert isinstance(j.out_cram, JobResourceFile)
            assert isinstance(j.out_crai, JobResourceFile)
            j.out_cram.add_extension('.cram')
            j.out_crai.add_extension('.cram.crai')

            b.write_output(j.out_cram, str(newpath))
            b.write_output(j.out_crai, str(newpath.with_suffix('.cram.crai')))

            print(f'Added reheader job for {path}')

        b.run(wait=False)

    else:
        sys.exit(__doc__)
