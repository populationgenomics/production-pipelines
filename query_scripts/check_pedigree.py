#!/usr/bin/env python3

"""
This script parses `somalier relate` (https://github.com/brentp/somalier) outputs,
and returns a non-zero code if either sex or pedigree mismatches the data in a provided PED file.
"""

import contextlib
import logging
import sys
import tempfile
from os.path import join, basename
from typing import Dict, List

import pandas as pd
import click
from peddy import Ped
from cpg_production_pipelines import utils

logger = logging.getLogger('check-pedigree')
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--somalier-samples',
    'somalier_samples_fpath',
    required=True,
    help='Path to somalier {prefix}.samples.tsv output file',
)
@click.option(
    '--somalier-pairs',
    'somalier_pairs_fpath',
    required=True,
    help='Path to somalier {prefix}.pairs.tsv output file',
)
@click.option(
    '--somalier-html',
    'somalier_html_fpath',
    help='Path to somalier {prefix}.html output file',
)
def main(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    somalier_html_fpath: str,
):  # pylint: disable=missing-function-docstring
    if somalier_samples_fpath.startswith('gs://'):
        local_tmp_dir = tempfile.mkdtemp()
        utils.gsutil_cp(somalier_samples_fpath, somalier_pairs_fpath)
        local_somalier_pairs_fpath = join(local_tmp_dir, basename(somalier_pairs_fpath))
        local_somalier_samples_fpath = join(
            local_tmp_dir, basename(somalier_samples_fpath)
        )
    else:
        local_somalier_pairs_fpath = somalier_pairs_fpath
        local_somalier_samples_fpath = somalier_samples_fpath

    logger.info('* Checking sex *')
    df = pd.read_csv(local_somalier_samples_fpath, delimiter='\t')
    missing_sex = df['original_pedigree_sex'] == 'unknown'
    mismatching_female = (df['sex'] == 2) & (df['original_pedigree_sex'] == 'male')
    mismatching_male = (df['sex'] == 1) & (df['original_pedigree_sex'] == 'female')
    mismatching_sex = mismatching_female | mismatching_male

    if (mismatching_sex | missing_sex).any():
        logger.info(f'Found PED samples with mismatching or missing sex:')
        for _, row in df[mismatching_sex | missing_sex].iterrows():
            inferred_sex = {1: 'male', 2: 'female'}.get(row.sex, 'unknown')
            logger.info(
                f'\t{row.sample_id} (provided: {row.original_pedigree_sex}, inferred: {inferred_sex})'
            )

    if not mismatching_sex.any():
        if missing_sex.any():
            logger.info(
                f'Inferred sex and pedigree matches or can be inferred for all samples.'
            )
        else:
            logger.info(f'Inferred sex and pedigree matches for all samples.')
    logger.info('-' * 10)

    logger.info('* Checking relatedness *')
    ped = Ped(local_somalier_samples_fpath)
    sample_by_id = {s.sample_id: s for s in ped.samples()}
    mismatching_pairs = []
    pairs_df = pd.read_csv(local_somalier_pairs_fpath, delimiter='\t')
    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        inferred_rel = infer_relationship(row['relatedness'], row['ibs0'], row['ibs2'])
        # Supressing all logging output from peddy as it would clutter the logs
        with contextlib.redirect_stderr(None):
            with contextlib.redirect_stdout(None):
                provided_rel = ped.relation(sample_by_id[s1], sample_by_id[s2])

        provided_to_inferred: Dict[str, List[str]] = {
            'unrelated': ['unrelated'],
            'related at unknown level': ['unrelated'],  # e.g. mom-dad
            'mom-dad': ['unrelated'],
            'parent-child': ['parent-child'],
            'grandchild': ['below_first_degree'],
            'niece/nephew': ['below_first_degree'],
            'great-grandchild': ['below_first_degree'],
            'cousins': ['below_first_degree'],
            'full siblings': ['siblings'],
            'siblings': ['siblings'],
            'unknown': [],
        }
        if inferred_rel not in provided_to_inferred[provided_rel]:
            mismatching_pairs.append(
                f'"{s1}" and "{s2}", '
                f'provided relationship: "{provided_rel}", inferred: "{inferred_rel}"'
            )
        pairs_df.loc[idx, 'provided_rel'] = provided_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    if mismatching_pairs:
        logger.info(f'Found sample pairs with mismatched relatedness:')
        for pair in mismatching_pairs:
            logger.info(f'\t{pair}')
    else:
        logger.info(f'Inferred pedigree matches for all samples.')

    logger.info('-' * 10)
    print_info(
        df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
        somalier_html_fpath,
    )
    if mismatching_sex.any() or mismatching_pairs:
        sys.exit(1)


def infer_relationship(coeff: float, ibs0: float, ibs2: float) -> str:
    """
    Inferres relashionship labels based on the kin coefficient
    and ibs0 and ibs2 values.
    """
    result = 'ambiguous'
    if coeff < 0.1:
        result = 'unrelated'
    elif 0.1 <= coeff < 0.38:
        result = 'below_first_degree'
    elif 0.38 <= coeff <= 0.62:
        if ibs0 / ibs2 < 0.005:
            result = 'parent-child'
        elif 0.015 < ibs0 / ibs2 < 0.052:
            result = 'siblings'
        else:
            result = 'first_degree'
    elif coeff > 0.8:
        result = 'duplicate_or_twins'
    return result


def print_info(
    samples_df,
    pairs_df,
    somalier_samples_fpath,
    somalier_pairs_fpath,
    somalier_html_fpath,
):
    """
    Print useful information to manually review pedigree check results
    """
    samples_str = samples_df.to_string()
    pairs_str = pairs_df[
        [
            '#sample_a',
            'sample_b',
            'relatedness',
            'ibs0',
            'ibs2',
            'n',
            'expected_relatedness',
        ]
    ].to_string()
    logger.info('')
    logger.info(
        f'Somalier results, samples (based on {somalier_samples_fpath}):\n{samples_str}\n'
    )
    logger.info(
        f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n{pairs_str}\n'
    )
    if somalier_html_fpath:
        logger.info(f'Somalier HTML report: {somalier_html_fpath}\n')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
