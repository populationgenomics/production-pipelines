#!/usr/bin/env python3

import hashlib
import os
import subprocess as sb
import sys

import click

from hailtop.batch.resource import JobResourceFile

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, get_config, image_path
from metamist import models
from metamist.apis import AnalysisApi
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


def query_metamist(dataset: str):
    cram_list = []

    for seqgroup in query(FIND_CRAMS, {'project': dataset})['project']['sequencingGroups']:
        for analysis in seqgroup['analyses']:
            fname = analysis['output']
            if fname is not None and fname.endswith('.cram'):
                cram_list.append(fname)
            else:
                sgid = seqgroup['id']
                anid = analysis['id']
                print(f'Sequencing group {sgid} analysis {anid} is not CRAM: {fname}')

    return cram_list


# Precomputed hashes (as per classify_SQ()) for common reference sets
KNOWN_REFS = {
    '49f7a7583b832aed0104b1901f9cf5fd': 'Homo_sapiens_assembly38_masked',
    '760a471d32c9a356ffb37f632bebfea3': 'Homo_sapiens_assembly38[unmasked]',
    }


def parse_one_SQ(line: str):
    table = {}
    for field in line.strip().split('\t')[1:]:
        key, value = field.split(':', maxsplit=1)
        table[key] = value
    return table


def classify_SQ(path: str):
    """Reads the file's headers and recognises well-known reference collections"""
    run = sb.run(['samtools', 'head', path], capture_output=True, text=True)
    if run.returncode != 0:
        return None  # Probably file not found

    m5 = {}
    for line in run.stdout.splitlines():
        if line.startswith('@SQ'):
            field = parse_one_SQ(line)
            m5[field['SN']] = field['M5']

    text = '/'.join([f'{name}:{m5[name]}' for name in sorted(m5.keys())])
    hashed = hashlib.md5(text.encode('ASCII')).hexdigest()

    return KNOWN_REFS.get(hashed, hashed)


def do_reheader(dry_run, refset, newref, newrefset, incram, outcram, outcrai):
    ref = classify_SQ(incram)
    if ref != refset:
        raise ValueError(f'{incram}: unexpected original reference headers: {ref}')

    translate = {}
    with open(newref) as refp:
        for line in refp:
            if line.startswith('@SQ'):
                field = parse_one_SQ(line)
                translate[field['SN']] = (line.strip(), field['LN'])

    head_cmd = ['samtools', 'head', incram]
    headers = sb.run(head_cmd, check=True, capture_output=True, text=True).stdout

    newhdr_file = os.path.join(os.environ['BATCH_TMPDIR'], 'tmp.hdr')
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

    if dry_run:
        print(f'{incram} original headers: {ref}')
        print(f'Would rewrite it to {outcram} with headers:')
        with open(newhdr_file) as fp:
            print(fp.read())
        return {}

    sb.run(['samtools', 'reheader', '--no-PG', '--in-place', newhdr_file, incram], check=True)
    os.rename(incram, outcram)
    sb.run(['samtools', 'index', '-o', outcrai, outcram], check=True)

    newref = classify_SQ(outcram)
    print(f'{outcram}: reheadered: {newref}')
    if newrefset is not None and newref != newrefset:
        raise ValueError(f'{outcram}: unexpected reference headers written: {newref}')

    return {
        'cram_size': os.stat(outcram).st_size,
    }


def do_metamist_update(dataset, sequencing_type, seqgroup, oldpath, newpath, result):
    analysis = models.Analysis(
        type='cram',
        status=models.AnalysisStatus('completed'),
        output=str(newpath),
        sequencing_group_ids=[seqgroup],
        meta={
            'sequencing_type': sequencing_type,
            'size': result['cram_size'],
            'source': f'Reheadered from {oldpath}',
        },
    )

    aid = AnalysisApi().create_analysis(dataset, analysis)
    print(f'Created Analysis(id={aid}, output={newpath}) in {dataset}')


@click.command()
@click.argument('source', nargs=-1)
@click.option('--ref', 'new_reference_path', required=True, metavar='GCS_PATH', help='Dict file containing the new @SQ headers')
@click.option('--expected', 'expected_refset', required=True, metavar='NAME', help='KNOWN_REFS name of the (bad) refset expected')
@click.option('--intended', 'new_refset', metavar='NAME', help='KNOWN_REFS name of the replacement headers')
@click.option('--dry-run', is_flag=True, help='Display information only without making changes')
def main(
    source: tuple[str],
    new_reference_path: str,
    expected_refset: str,
    new_refset: str,
    dry_run: bool,
):
    """
    This script uses the usual config entries (workflow.dataset, .access_level, .sequencing_type
    and images.samtools) and the following options to select CRAM files to have their @SQ
    headers rewritten.
    """
    config = get_config(True)

    dataset = config['workflow']['dataset']
    if config['workflow']['access_level'] == 'test' and not dataset.endswith('-test'):
        dataset = f'{dataset}-test'

    if len(source) > 0:
        print(f'Scanning {" ".join(source)} for *.cram files')
        cram_list = [fn for srcdir in source
                     for fn in to_path(srcdir).iterdir() if fn.suffix == '.cram']
    else:
        print(f'Finding CRAMs for {dataset} in database')
        cram_list = query_metamist(dataset)

    print(f'Total CRAM file entries found: {len(cram_list)}')

    sequencing_type = config['workflow']['sequencing_type']

    reheader_list = []
    for fname in cram_list:
        path = to_path(fname)
        reheadered_path = path.parent / 'reheadered' / path.name

        if not path.exists():
            print(f'{fname}: file does not exist')
        elif reheadered_path.exists():
            print(f'{reheadered_path}: reheadered file already exists')
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
        j.storage(filesize * 1.2)  # Allow some extra space for index file, references, etc

        j_result = j.call(do_reheader,
                          dry_run, expected_refset, new_reference, new_refset,
                          b.read_input(str(path)), j.out_cram, j.out_crai)

        assert isinstance(j.out_cram, JobResourceFile)
        assert isinstance(j.out_crai, JobResourceFile)
        j.out_cram.add_extension('.cram')
        j.out_crai.add_extension('.cram.crai')

        if not dry_run:
            b.write_output(j.out_cram, str(newpath))
            b.write_output(j.out_crai, str(newpath.with_suffix('.cram.crai')))

        print(f'Added reheader job for {path}')

        if not dry_run:
            db_j = b.new_python_job(f'Update metamist for {newpath}')
            db_j.image(config['workflow']['driver_image'])

            db_j.call(do_metamist_update,
                      dataset, sequencing_type, path.stem, str(path), str(newpath), j_result)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter  # click will add the arguments
