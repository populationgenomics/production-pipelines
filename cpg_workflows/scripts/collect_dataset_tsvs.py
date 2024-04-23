#!/usr/bin/env python3

"""
Collect all the dataset TSVs and save them in a single file
"""

import re
from argparse import ArgumentParser
from csv import DictReader

# detect "OMIM:12345" or "ORPHA:12345"
DISEASE_RE = re.compile(r'((?:OMIM|ORPHA):\d+)')

# at time of writing
FIELDNAMES = [
    '#RANK',
    'ID',
    'GENE_SYMBOL',
    'ENTREZ_GENE_ID',
    'MOI',
    'P-VALUE',
    'EXOMISER_GENE_COMBINED_SCORE',
    'EXOMISER_GENE_PHENO_SCORE',
    'EXOMISER_GENE_VARIANT_SCORE',
    'HUMAN_PHENO_SCORE',
    'MOUSE_PHENO_SCORE',
    'FISH_PHENO_SCORE',
    'WALKER_SCORE',
    'PHIVE_ALL_SPECIES_SCORE',
    'OMIM_SCORE',
    'MATCHES_CANDIDATE_GENE',
    'HUMAN_PHENO_EVIDENCE',
    'MOUSE_PHENO_EVIDENCE',
    'FISH_PHENO_EVIDENCE',
    'HUMAN_PPI_EVIDENCE',
    'MOUSE_PPI_EVIDENCE',
    'FISH_PPI_EVIDENCE',
]


def _parse_tsv_row(row):
    """stolen from seqr"""
    return [s.strip('#').strip().strip('"') for s in row.rstrip('\n').split('\t')]


def main(project: str, input_files: list[str], p_threshold: float = 0.05):
    """
    Collect all the dataset TSVs and save them in a single file
    Embellish the data with the mapping and project name

    Args:
        project ():
        input_files (): all input files to process
        p_threshold (): above this value, ignore the row and break
    """

    # these are the headings for Seqr
    output_lines: list[list[str]] = [
        [
            'tool',
            'project',
            'sampleId',
            'rank',
            'geneId',
            'diseaseId',
            'diseaseName',
            'scoreName1',
            'score1',
            'scoreName2',
            'score2',
            'scoreName3',
            'score3',
        ],
    ]

    for family_file in input_files:
        print(f'Processing {family_file}')

        # crack it open, and take a sip
        with open(family_file, 'r') as input_handle:
            reader = DictReader(input_handle, delimiter='\t')
            for row in reader:
                if float(row['P-VALUE']) > p_threshold:
                    break

                # fish out the relevant disease details, if present
                disease_name: str = row['HUMAN_PHENO_EVIDENCE'].split(' (')[0]
                disease_id = ''
                if result := re.search(DISEASE_RE, row['HUMAN_PHENO_EVIDENCE']):
                    disease_id = result.group(0)

                # get participant_id from the filename
                participant_id = family_file.split('/')[-1].split('.')[0]

                output_lines.append(
                    [
                        'exomiser',
                        project,
                        participant_id,
                        row['#RANK'],
                        row['GENE_SYMBOL'],
                        disease_id,
                        disease_name,
                        'exomiser_score',
                        row['EXOMISER_GENE_COMBINED_SCORE'],
                        'phenotype_score',
                        row['EXOMISER_GENE_PHENO_SCORE'],
                        'variant_score',
                        row['EXOMISER_GENE_VARIANT_SCORE'],
                    ],
                )
    return output_lines


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--project', help='Name of Seqr Project')
    parser.add_argument('--input', help='All input files', nargs='+')
    parser.add_argument('--output', help='Output file to save the combined data')
    parser.add_argument('--p_threshold', type=float, help='P-value threshold for filtering', default=0.05)
    args = parser.parse_args()

    data_lines = main(args.project, args.input, args.p_threshold)
    print(f'Found {len(data_lines)} lines of data')

    with open(args.output, 'w') as handle:
        for line in data_lines:
            handle.write('\t'.join(line) + '\n')
