#!/usr/bin/env python3

"""
Copy CRAM files while replacing existing @SQ headers with new ones as specified,
add new analysis records to metamist pointing to the new CRAM files, and
inactivate the existing analysis records for these sequencing groups.

Typical usage:

analysis-runner ... scripts/fix_sq_headers --ref gs://.../reference.dict
    --expected 'Homo_sapiens_assembly38[unmasked]'
    --intended 'Homo_sapiens_assembly38_masked'
    [gs://... directory or file paths to limit processing]...
"""

import hashlib
import os
import subprocess as sb
import sys

import click

from hailtop.batch.resource import JobResourceFile

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import get_batch
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
      sequencingGroups {
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
                    raise ValueError(
                        f'{incram}:{field["SN"]}: new length {new_len} differs from existing {field["LN"]}',
                    )
                print(new_line, file=newhdr)
            else:
                print(line, file=newhdr)

    if dry_run:
        print(f'{incram} original headers: {ref}')
        print(f'Would rewrite it to {outcram} with headers:')
        with open(newhdr_file) as fp:
            print(fp.read())
        return {}

    sb.run(
        ['samtools', 'reheader', '--no-PG', '--in-place', newhdr_file, incram],
        check=True,
    )
    os.rename(incram, outcram)
    sb.run(['samtools', 'index', '-o', outcrai, outcram], check=True)

    newref = classify_SQ(outcram)
    print(f'{outcram}: reheadered: {newref}')
    if newrefset is not None and newref != newrefset:
        raise ValueError(f'{outcram}: unexpected reference headers written: {newref}')

    return {
        'cram_size': os.stat(outcram).st_size,
    }


def do_metamist_update(dataset, activeseqgroup, oldanalysis, oldpath, newpath, result):
    aapi = AnalysisApi()

    newmeta = oldanalysis['meta']
    newmeta['size'] = result['cram_size']
    newmeta['source'] += f' (reheadered from {oldpath})'

    newanalysis = models.Analysis(
        type=oldanalysis['type'],
        status=models.AnalysisStatus(oldanalysis['status'].lower()),
        output=str(newpath),
        sequencing_group_ids=[activeseqgroup['id']],
        meta=newmeta,
    )
    aid = aapi.create_analysis(dataset, newanalysis)
    print(f'Created Analysis(id={aid}, output={newpath}) in {dataset}')

    for analysis in activeseqgroup['analyses']:
        activeid = analysis['id']
        analysisupdate = models.AnalysisUpdateModel(
            active=False,
            status=models.AnalysisStatus(analysis['status'].lower()),
        )
        if not aapi.update_analysis(activeid, analysisupdate):
            raise ValueError(f'Updating analysis id={activeid} failed')
        print(f'Inactivated Analysis(id={activeid}, output={analysis["output"]} in {dataset}')


@click.command()
@click.argument('source', nargs=-1)
@click.option(
    '--ref',
    'new_reference_path',
    required=True,
    metavar='GCS_PATH',
    help='Dict file containing the new @SQ headers',
)
@click.option(
    '--expected',
    'expected_refset',
    required=True,
    metavar='NAME',
    help='KNOWN_REFS name of the (bad) refset expected',
)
@click.option(
    '--intended',
    'new_refset',
    metavar='NAME',
    help='KNOWN_REFS name of the replacement headers',
)
@click.option('--dry-run', is_flag=True, help='Display information only without making changes')
def main(
    source: tuple[str],
    new_reference_path: str,
    expected_refset: str,
    new_refset: str,
    dry_run: bool,
):
    """
    This script uses the usual config entries (workflow.dataset, .access_level, .driver_image
    and images.samtools) and the following options to select CRAM files to have their @SQ
    headers rewritten.
    """
    config = get_config(True)

    dataset = config['workflow']['dataset']
    if config['workflow']['access_level'] == 'test' and not dataset.endswith('-test'):
        dataset = f'{dataset}-test'

    project_analyses = query(FIND_ANALYSES, {'project': dataset})['project']['analyses']
    analysis_by_path = {analysis['output']: analysis for analysis in project_analyses}

    if len(source) > 0:
        print(f'Scanning {" ".join(source)} for *.cram files')
        cram_list = []
        for src in source:
            if src.endswith('.cram'):
                if src in analysis_by_path:
                    cram_list.append(src)
                else:
                    print(f'No analysis record found for {src}')
            else:
                for fn in to_path(src).iterdir():
                    if fn.suffix == '.cram':
                        if fn in analysis_by_path:
                            cram_list.append(fn)
                        else:
                            print(f'No analysis record found for {fn}')
    else:
        print(f'Using all CRAMs for {dataset} in database')
        cram_list = analysis_by_path.keys()

    print(f'Total CRAM file entries found: {len(cram_list)}')

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

    project_samples = query(FIND_SAMPLES, {'project': dataset})['project']['samples']
    sample_by_id = {sample['id']: sample for sample in project_samples}

    b = get_batch(f'Reheader {dataset} to {new_reference_path}')
    new_reference = b.read_input(new_reference_path)

    for path, newpath, filesize in reheader_list:
        j = b.new_python_job(f'Reheader {path}')
        j.image(image_path('samtools'))
        j.storage(filesize * 1.2)  # Allow some extra space for index file, references, etc

        j_result = j.call(
            do_reheader,
            dry_run,
            expected_refset,
            new_reference,
            new_refset,
            b.read_input(str(path)),
            j.out_cram,
            j.out_crai,
        )

        assert isinstance(j.out_cram, JobResourceFile)
        assert isinstance(j.out_crai, JobResourceFile)
        j.out_cram.add_extension('.cram')
        j.out_crai.add_extension('.cram.crai')

        if not dry_run:
            b.write_output(j.out_cram, str(newpath))
            b.write_output(j.out_crai, str(newpath.with_suffix('.cram.crai')))

        print(f'Added reheader job for {path}')

        analysis = analysis_by_path[str(path)]
        if seqgroupcount := len(analysis['sequencingGroups']) != 1:
            print(f'Analysis id={analysis["id"]} for {path} is in {seqgroupcount} sequencing groups')
        elif not dry_run:
            db_j = b.new_python_job(f'Update metamist for {newpath}')
            db_j.image(config['workflow']['driver_image'])

            sample = sample_by_id[analysis['sequencingGroups'][0]['sample']['id']]
            seqgroup = sample['sequencingGroups'][0]

            db_j.call(
                do_metamist_update,
                dataset,
                seqgroup,
                analysis,
                str(path),
                str(newpath),
                j_result,
            )

            oldids = ', '.join(str(a['id']) for a in seqgroup['analyses'])
            print(f'Added metamist job to create new analysis and inactivate id={oldids}')

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter  # click will add the arguments
