#!/usr/bin/env python3

"""
Change the reference contigs md5-sums in a CRAM header 
from alts-masked to alts-unmasked versions.

Usage:

```sh
samtools reheader sample.cram --in-place \
    --command "fix_cram_header.py Homo_sapiens_assembly38.dict"
```
"""
import sys


def _parse_sq_line(line: str) -> dict[str, str] | None:
    if not line.startswith('@SQ'):
        return None
    cols = line.strip().split('\t')[1:]
    return {(p := col.split(':', maxsplit=1))[0]: p[1] for col in cols}


def main():
    fasta_dict_path = sys.argv[1]
    unmasked_dict = dict()
    with open(fasta_dict_path) as f:
        for line in f:
            if record := _parse_sq_line(line):
                unmasked_dict[record['SN']] = record

    for line in sys.stdin:
        line = line.strip()
        if record := _parse_sq_line(line):
            unmasked_record = unmasked_dict[record['SN']]
            if unmasked_record['LN'] != record['LN']:
                print(f'{record["SN"]} lengths do not match:', file=sys.stderr)
                print(f'original record: {record}', file=sys.stderr)
                print(f'unmasked record: {unmasked_record}', file=sys.stderr)

            record['M5'] = unmasked_record['M5']  # fixing MD5 sums
            line = '\t'.join(['@SQ'] + [f'{k}:{v}' for k, v in record.items()])
            print(line)
        else:
            print(line)


if __name__ == '__main__':
    main()
