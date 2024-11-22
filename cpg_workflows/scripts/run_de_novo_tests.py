"""
do some de novo testing, to see how many of our knowns we detect

analysis-runner --dataset seqr \
  --description 'run de novo testing' \
  --access-level test \
  --output-dir denovo_testing \
  python3 cpg_workflows/scripts/run_de_novo_tests.py \
    --mt gs://cpg-seqr-test/de_novo_experimentation/repartitioned_filtered_de_novos.mt \
    --ped gs://cpg-seqr-test/de_novo_experimentation/de_novo_fams.ped \
    --out gs://cpg-seqr-test/de_novo_experimentation/results
"""

import logging
from argparse import ArgumentParser
from itertools import product

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch


def main(mt_path: str, pedigree: str, destination: str):
    """
    run de novo tests using multiple parameter combinations

    Args:
        mt_path ():
        pedigree ():
        destination ():
    """

    # read in the MT of isolated de novo variant calls, and the pedigree
    mt = hl.read_matrix_table(mt_path)
    pedigree = hl.Pedigree.read(pedigree)
    destination_path = to_path(destination)

    # upgrade PL to be perfect if missing - this is a major factor in success rate
    # we do this in Talos due to multiple misses in the past
    if correct_pls := config_retrieve(['de_novo', 'correct_pls'], True):
        mt = mt.annotate_entries(
            PL=hl.case()
            .when(~hl.is_missing(mt.PL), mt.PL)
            .when((mt.GT.is_non_ref()) | (hl.is_missing(mt.GQ)), hl.missing('array<int32>'))
            .default([0, mt.GQ, 1000]),
        )

    # pull these variables from config - default to the real hail defaults
    # to expand this testing we can add elements to the list
    # this script will test every combination of these parameters
    min_gqs = config_retrieve(['de_novo', 'min_gqs'], [20])
    min_ps = config_retrieve(['de_novo', 'min_ps'], [0.05])
    max_parent_abs = config_retrieve(['de_novo', 'max_parent_abs'], [0.05])
    min_child_abs = config_retrieve(['de_novo', 'min_child_abs'], [0.2])
    min_dp_ratios = config_retrieve(['de_novo', 'min_dp_ratios'], [0.1])

    for min_gq, min_p, max_parent_ab, min_child_ab, min_dp_ratio in product(
        min_gqs,
        min_ps,
        max_parent_abs,
        min_child_abs,
        min_dp_ratios,
    ):
        logging.info(min_gq, min_p, max_parent_ab, min_child_ab, min_dp_ratio, correct_pls)
        args_as_str = '_'.join([str(x) for x in [min_gq, min_p, max_parent_ab, min_child_ab, min_dp_ratio]])
        dest_path = destination_path / f'{args_as_str}.ht'
        if (dest_path / '_SUCCESS').exists():
            logging.info(f'Skipping {args_as_str} as it has already been run')
            continue

        dn_table = hl.de_novo(
            mt,
            pedigree,
            pop_frequency_prior=mt.gnomad_exomes.AF,
            min_gq=min_gq,
            min_p=min_p,
            max_parent_ab=max_parent_ab,
            min_child_ab=min_child_ab,
            min_dp_ratio=min_dp_ratio,
            ignore_in_sample_allele_frequency=True,
        )

        successes = dn_table.count()
        logging.info(f'Successes: {successes}')
        if successes < 5:
            logging.info('Too few successes, skipping write')
            continue

        # write the table out
        dn_table.write(str(dest_path))
        # then write an additional file with the number of wins
        with open(str(dest_path / 'num_wins.txt'), 'w') as f:
            f.write(str(successes))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--mt', help='MatrixTable')
    parser.add_argument('--ped', help='Pedigree')
    parser.add_argument('--out', help='Destination to write output to')
    args = parser.parse_args()

    # boot up Hail
    init_batch()

    # turn up the standard logging
    logging.basicConfig(level=logging.INFO)

    main(
        mt_path=args.mt,
        pedigree=args.ped,
        destination=args.out,
    )
