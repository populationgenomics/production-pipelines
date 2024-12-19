"""
script for taking the per-family TSVs from Exomiser and combining them
into a single JSON file, also written as a Hail Table
"""

import json
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader

import hail as hl

from cpg_utils import to_path

ORDERED_ALLELES: list[str] = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']
FAM_DICT = dict[str, dict[str, list[dict]]]
VAR_DICT = dict[str, list[str]]


def process_tsv(tsv_path: str) -> tuple[FAM_DICT, VAR_DICT]:
    """
    Read a single TSV file, and return the parsed results
    Return is in two forms:
    - a dictionary of variant IDs to lists of sample IDs & details


    the return is a dictionary of variant keys to lists of details
    the same variant can apply to the same family with multiple MOIs so we retain all

    Args:
        tsv_path (str): path to this TSV file

    Returns:
        a dictionary of variant keys to lists of the 'family_rank_moi'
    """

    # this is on the assumption that we retain the path/to/file/FAMILY.variant.tsv naming convention
    # if we read into the batch as /tmp/path/FAMILY this logic is still valid
    file_as_path = to_path(tsv_path)
    family = file_as_path.name.split('.')[0]

    variant_dictionary: VAR_DICT = defaultdict(list)
    family_dictionary: FAM_DICT = {family: defaultdict(list)}

    with file_as_path.open() as handle:
        reader = DictReader(handle, delimiter='\t')
        for row in reader:
            assert isinstance(row, dict)

            # this isn't prefixed with 'chr' in the TSV, so we add it
            chrom = f'chr{row["CONTIG"]}'
            pos = row['START']
            ref = row['REF']
            alt = row['ALT']

            variant_key = f'{chrom}:{pos}:{ref}:{alt}'

            rank = int(row['#RANK'])
            moi = row['MOI']

            # adds easily parsed details to the family dictionary, easily extensible
            family_dictionary[family][variant_key].append(
                {
                    'rank': rank,
                    'moi': moi,
                },
            )

            # compressed representation for exporting to Hail
            squashed_details = f'{family}_{rank}_{moi}'
            variant_dictionary[variant_key].append(squashed_details)

    return family_dictionary, variant_dictionary


def process_and_sort_variants(all_variants: VAR_DICT) -> list[dict]:
    """
    applies dual-layer sorting to the list of all variants

    Args:
        all_variants (): list of all submissions

    Returns:
        a list of variants, sorted hierarchically on chr & pos
    """

    # first, split the data up into component parts
    all_vars: list[dict] = []
    for variant_key, family_details in all_variants.items():
        chrom, pos, ref, alt = variant_key.split(':')
        all_vars.append(
            {'contig': chrom, 'position': int(pos), 'alleles': [ref, alt], 'family_details': '::'.join(family_details)},
        )

    # then sort on chr and position
    return sorted(all_vars, key=lambda x: (ORDERED_ALLELES.index(x['contig']), x['position']))


def munge_into_hail_table(all_variants: VAR_DICT, output_path: str):
    """
    Reformat and sort the data into a Hail Table

    Args:
        all_variants ():
        output_path ():
    """

    # process the data into a list of dicts
    var_data = process_and_sort_variants(all_variants)

    # write it out to a local temp file
    temp_filepath = 'temporary_out.json'
    with open(temp_filepath, 'w', encoding='utf-8') as handle:
        for each_dict in var_data:
            handle.write(f'{json.dumps(each_dict)}\n')

    # start a hail runtime, local style
    hl.context.init_local(default_reference='GRCh38')

    # define the schema for each written line
    schema = hl.dtype(
        'struct{' 'contig:str,' 'position:int32,' 'alleles:array<str>,' 'family_details:str' '}',
    )

    # import the table, and transmute to top-level attributes
    ht = hl.import_table(temp_filepath, no_header=True, types={'f0': schema})
    ht = ht.transmute(**ht.f0)

    # create a locus value, and key the table by this
    ht = ht.transmute(locus=hl.locus(ht.contig, ht.position))
    ht = ht.key_by(ht.locus, ht.alleles)

    # write out to the specified location
    ht.write(f'{output_path}.ht', overwrite=True)


def main(input_tsvs: list[str], output_path: str, as_hail: bool = True):
    """
    Combine the per-family TSVs into a single JSON file, and write as a Hail Table

    The data will be aggregated per variant, not per family

    Args:
        input_tsvs (list[str]): list of paths to the per-family TSVs
        output_path (str): where to write the output
        as_hail (bool): if True, write the data out as a Hail Table
    """

    variant_dictionary: VAR_DICT = defaultdict(list)
    all_family_dictionary: FAM_DICT = {}
    for tsv in input_tsvs:
        family_dict, variant_dict = process_tsv(tsv)
        all_family_dictionary.update(family_dict)
        for key, value_list in variant_dict.items():
            variant_dictionary[key].extend(value_list)

    # write the JSON
    with open(f'{output_path}.json', 'w') as handle:
        json.dump(all_family_dictionary, handle, indent=4)

    if as_hail:
        munge_into_hail_table(variant_dictionary, output_path)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the per-family TSVs', nargs='+', required=True)
    parser.add_argument('--output', help='Where to write the output, extended as .json and .ht', required=True)
    args = parser.parse_args()

    main(input_tsvs=args.input, output_path=args.output)


if __name__ == '__main__':
    cli_main()
