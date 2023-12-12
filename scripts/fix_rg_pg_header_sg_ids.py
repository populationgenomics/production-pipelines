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
query CramQuery($project: String!, $sequencing_groups: [String!]) {
  project(name: $project) {
    sequencingGroups(id: {in_:$sequencing_groups}) {
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


def query_metamist(dataset: str, sequencing_groups: list[str]) -> dict[str, str]:
    cram_list = {}

    for seqgroup in query(FIND_CRAMS, {'project': dataset, 'sequencing_groups': sequencing_groups})['project']['sequencingGroups']:
        sgid = seqgroup['id']
        for analysis in seqgroup['analyses']:
            fname = analysis['output']
            if fname is not None and fname.endswith('.cram'):
                cram_list[sgid] = fname
            else:
                anid = analysis['id']
                print(f'Sequencing group {sgid} analysis {anid} is not CRAM: {fname}')

    return cram_list


def parse_one_header_line(line: str):
    table = {}
    for field in line.strip().split('\t')[1:]:
        key, value = field.split(':', maxsplit=1)
        table[key] = value
    return table


def extract_header_sgid(path: str):
    """Reads the file's headers and finds the sequencing group IDs within"""
    run = sb.run(['samtools', 'head', path], capture_output=True, text=True)
    if run.returncode != 0:
        return None  # Probably file not found

    header_sgid = {}
    for line in run.stdout.splitlines():
        if line.startswith('@RG'):
            rg_field = parse_one_header_line(line)
            header_sgid['ID'] = rg_field['ID']
            header_sgid['SM'] = rg_field['SM']
            assert rg_field['ID'] == rg_field['SM']
        if line.startswith('@PG'):
            pg_field = parse_one_header_line(line)
            cmd = pg_field['CL']
            if cmd[:-4] == ['--RGID', header_sgid['ID'], '--RGSM', header_sgid['ID']]:
                header_sgid['RGID'] = cmd[-3]
                header_sgid['RGSM'] = cmd[-1]
                assert header_sgid['RGID'] == header_sgid['RGSM']
                assert header_sgid['RGID'] == header_sgid['ID']

    return header_sgid


def do_reheader(dry_run, sgid, incram, outcram, outcrai):
    header = extract_header_sgid(incram)
    # Check if any of the values in the header dict are not the sgid
    if not any(header.get(key) != sgid for key in header):
        raise ValueError(f'{incram}: headers already match {sgid}')

    head_cmd = ['samtools', 'head', incram]
    headers = sb.run(head_cmd, check=True, capture_output=True, text=True).stdout

    newhdr_file = os.path.join(os.environ['BATCH_TMPDIR'], 'tmp.hdr')
    with open(newhdr_file, 'w') as newhdr:
        for line in headers.splitlines():
            pg_line_count = 0
            if line.startswith('@RG'):
                new_line = line.replace(header['ID'], sgid)
                print(new_line, file=newhdr)
            elif line.startswith('@PG') and pg_line_count < 1:
                new_line = line.replace(header['ID'], sgid)
                print(line, file=newhdr)
                pg_line_count += 1      
            else:
                print(line, file=newhdr)

    if dry_run:
        print(f'{incram} original header sequencing group ID: {header["ID"]}')
        print(f'Would rewrite it to {outcram} with headers:')
        with open(newhdr_file) as fp:
            print(fp.read())
        return {}

    sb.run(['samtools', 'reheader', '--no-PG', '--in-place', newhdr_file, incram], check=True)
    os.rename(incram, outcram)
    sb.run(['samtools', 'index', '-o', outcrai, outcram], check=True)

    new_header_sgid = extract_header_sgid(outcram)
    print(f'{outcram}: reheadered: {new_header_sgid}')

    return {
        'cram_size': os.stat(outcram).st_size,
    }


def do_metamist_update(dataset, sequencing_type, seqgroup, oldpath, newpath, result):
    aapi = AnalysisApi()

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
    aid = aapi.create_analysis(dataset, analysis)
    print(f'Created Analysis(id={aid}, output={newpath}) in {dataset}')


@click.command()
@click.option('--sequencing-groups', multiple=True, help='Sequencing group ID and name')
@click.option('--dry-run', is_flag=True, help='Display information only without making changes')
def main(
    sequencing_groups: list[str],
    dry_run: bool,
):
    """
    This script uses the usual config entries (workflow.dataset, .access_level, .sequencing_type
    and images.samtools) to check the CRAM headers for a list of sequencing groups, reheadering and
    renaming the files with the correct sequencing group.
    """
    config = get_config(True)

    dataset = config['workflow']['dataset']
    if config['workflow']['access_level'] == 'test' and not dataset.endswith('-test'):
        dataset = f'{dataset}-test'

    print(f'Finding CRAMs for {len(sequencing_groups)} sequencing groups from {dataset} in database')
    cram_list = query_metamist(dataset, sequencing_groups)

    print(f'Total CRAM file entries found: {len(cram_list)}')

    sequencing_type = config['workflow']['sequencing_type']

    reheader_list = []
    for sgid, fname in cram_list.items():
        if fname.endswith(f'{sgid}.cram'):
            continue
        
        path = to_path(fname)
        reheadered_path = path.parent / f'{sgid}.cram'

        if not path.exists():
            print(f'{fname}: file does not exist')
        elif reheadered_path.exists():
            print(f'{reheadered_path}: reheadered file already exists')
        else:
            reheader_list.append((path, reheadered_path, path.stat().st_size))

    if len(reheader_list) == 0:
        sys.exit('No CRAMs to be reheadered!')

    print(f'CRAM files to be reheadered: {len(reheader_list)}')

    b = get_batch(f'Reheader CRAMs in {dataset}')

    for path, newpath, filesize in reheader_list:
        j = b.new_python_job(f'Reheader {path}')
        j.image(image_path('samtools'))
        j.storage(filesize * 1.2)  # Allow some extra space for index file, references, etc

        j_result = j.call(do_reheader,
                          dry_run, reheadered_path.stem,
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
                      dataset, sequencing_type, reheadered_path.stem, str(path), str(newpath), j_result)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter  # click will add the arguments
