#!/usr/bin/env python3

"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a non-zero code if either sex or pedigree mismatches the data in a 
provided PED file.
"""

import contextlib
import logging
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from os.path import join, basename
from typing import Optional, List

import pandas as pd
import click
from peddy import Ped

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
    somalier_html_fpath: Optional[str],
):  # pylint: disable=missing-function-docstring
    if somalier_samples_fpath.startswith('gs://'):
        local_tmp_dir = tempfile.mkdtemp()
        cmd = f'gsutil cp {somalier_samples_fpath} {somalier_pairs_fpath}'
        subprocess.run(cmd, check=False, shell=True)
        local_somalier_pairs_fpath = join(local_tmp_dir, basename(somalier_pairs_fpath))
        local_somalier_samples_fpath = join(
            local_tmp_dir, basename(somalier_samples_fpath)
        )
    else:
        local_somalier_pairs_fpath = somalier_pairs_fpath
        local_somalier_samples_fpath = somalier_samples_fpath

    df = pd.read_csv(local_somalier_samples_fpath, delimiter='\t')
    df.sex = df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    df.original_pedigree_sex = df.original_pedigree_sex.apply(
        lambda x: {'-9': 'unknown'}.get(x, x)
    )
    bad_samples = list(df[df.gt_depth_mean == 0.0].sample_id)
    if bad_samples:
        logger.info(
            f'Excluding samples with non enough coverage to make inference: '
            f'{", ".join(bad_samples)}'
        )
    logger.info('-' * 10)

    logger.info('* Checking sex *')
    missing_inferred_sex = df.sex == 'unknown'
    missing_provided_sex = df.original_pedigree_sex == 'unknown'
    mismatching_female = (df.sex == 'female') & (df.original_pedigree_sex == 'male')
    mismatching_male = (df.sex == 'male') & (df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male

    def _print_stats(df_filter):
        for _, row in df[df_filter].iterrows():
            logger.info(
                f'\t{row.sample_id} ('
                f'provided: {row.original_pedigree_sex}, '
                f'inferred: {row.sex}, '
                f'mean depth: {row.gt_depth_mean})'
            )

    if mismatching_sex.any():
        logger.info(f'Found PED samples with mismatching sex:')
        _print_stats(mismatching_sex)
    if missing_inferred_sex.any():
        logger.info(f'Samples with missing inferred sex:')
        _print_stats(missing_inferred_sex)
    if missing_provided_sex.any():
        logger.info(f'Samples with missing provided sex:')
        _print_stats(missing_provided_sex)

    if not mismatching_sex.any():
        if (missing_inferred_sex | missing_provided_sex).any():
            logger.info(
                f'Inferred sex and pedigree matches or can be inferred for all samples.'
            )
        else:
            logger.info(f'Inferred sex and pedigree matches for all samples.')
    logger.info('-' * 10)

    logger.info('* Checking relatedness *')
    ped = Ped(local_somalier_samples_fpath)
    sample_by_id = {s.sample_id: s for s in ped.samples()}
    
    pairs_provided_as_unrelated_but_inferred_related = []
    other_mismatching_pairs = []
    pairs_df = pd.read_csv(local_somalier_pairs_fpath, delimiter='\t')
    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_samples or s2 in bad_samples:
            continue

        inferred_rel = infer_relationship(row['relatedness'], row['ibs0'], row['ibs2'])
        # Supressing all logging output from peddy as it would clutter the logs
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            peddy_rel = ped.relation(sample_by_id[s1], sample_by_id[s2])

        def _match_peddy_with_inferred(peddy_rel):
            return {
                'unrelated': 'unrelated',
                'related at unknown level': 'unrelated',
                'mom-dad': 'unrelated',
                'parent-child': 'parent-child',
                'grandchild': 'below_first_degree',
                'niece/nephew': 'below_first_degree',
                'great-grandchild': 'below_first_degree',
                'cousins': 'below_first_degree',
                'full siblings': 'siblings',
                'siblings': 'siblings',
                'unknown': 'unknown',
            }.get(peddy_rel)
        

        if (
            _match_peddy_with_inferred(peddy_rel) == 'unknown'
            and inferred_rel != 'unknown'
            or _match_peddy_with_inferred(peddy_rel) == 'unrelated'
            and inferred_rel != 'unrelated'
        ):
            pairs_provided_as_unrelated_but_inferred_related.append(
                f'"{s1}" and "{s2}", '
                f'provided: "{peddy_rel}", '
                f'inferred: "{inferred_rel}" '
                f'(rel={row["relatedness"]})'
            )
        elif inferred_rel != _match_peddy_with_inferred(peddy_rel):
            other_mismatching_pairs.append(
                f'"{s1}" and "{s2}", '
                f'provided: "{peddy_rel}", '
                f'inferred: "{inferred_rel}" '
                f'(rel={row["relatedness"]})'
            )
        pairs_df.loc[idx, 'provided_rel'] = peddy_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    if pairs_provided_as_unrelated_but_inferred_related:
        logger.info(
            f'Found sample pairs that are provided as unrelated, but '
            f'inferred as related:'
        )
        for pair in pairs_provided_as_unrelated_but_inferred_related:
            logger.info(f'\t{pair}')
    if other_mismatching_pairs:
        logger.info(f'Found sample pairs with mismatched relatedness:')
        for pair in other_mismatching_pairs:
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
    if mismatching_sex.any() or other_mismatching_pairs:
        sys.exit(1)


def infer_relationship(coeff: float, ibs0: float, ibs2: float) -> str:
    """
    Inferres relashionship labels based on the kin coefficient
    and ibs0 and ibs2 values.
    """
    if coeff < 0.2:
        result = 'unrelated'
    elif coeff < 0.38:
        result = 'below_first_degree'
    elif coeff <= 0.62:
        if ibs0 / ibs2 < 0.005:
            result = 'parent-child'
        elif 0.015 < ibs0 / ibs2 < 0.052:
            result = 'siblings'
        else:
            result = 'first_degree'
    elif coeff < 0.8:
        result = 'first_degree_or_duplicate_or_twins'
    elif coeff >= 0.8:
        result = 'duplicate_or_twins'
    else:
        result = 'nan'
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
    if len(samples_df) < 400:
        samples_str = samples_df.to_string()
        logger.info('')
        logger.info(
            f'Somalier results, samples (based on {somalier_samples_fpath}):\n{samples_str}\n'
        )
    if len(pairs_df) < 400:
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
        logger.info(
            f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n{pairs_str}\n'
        )
    if somalier_html_fpath:
        logger.info(f'Somalier HTML report: {somalier_html_fpath}\n')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
