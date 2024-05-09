#!/usr/bin/env python3

"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""

import contextlib
import logging
import os
from typing import Optional

import click
import pandas as pd
from peddy import Ped

from cpg_utils import to_path
from cpg_utils.slack import send_message

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


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
    '--expected-ped',
    'expected_ped_fpath',
    required=True,
    help='Path to PED file with expected pedigree',
)
@click.option('--html-url', 'html_url', help='Somalier HTML URL')
@click.option('--dataset', 'dataset', help='Dataset name')
@click.option('--title', 'title', help='Report title')
@click.option(
    '--send-to-slack/--no-send-to-slack',
    'send_to_slack',
    help='Send log to Slack message, according to environment variables SLACK_CHANNEL and SLACK_TOKEN',
)
def main(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    html_url: Optional[str] = None,
    dataset: Optional[str] = None,
    title: Optional[str] = None,
    send_to_slack: bool = True,
):
    """
    Report pedigree inconsistencies, given somalier outputs.
    """
    run(
        somalier_samples_fpath=somalier_samples_fpath,
        somalier_pairs_fpath=somalier_pairs_fpath,
        expected_ped_fpath=expected_ped_fpath,
        html_url=html_url,
        dataset=dataset,
        title=title,
        send_to_slack=send_to_slack,
    )


_messages: list[str] = []


def info(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logging.info(msg)


def warning(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logging.info(msg)


def error(msg):
    """
    Record and forward.
    """
    _messages.append(msg)
    logging.error(msg)


def run(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    html_url: Optional[str] = None,
    dataset: Optional[str] = None,
    title: Optional[str] = None,
    send_to_slack: bool = True,
):
    """
    Report pedigree inconsistencies, given somalier outputs.
    """
    logging.info(os.getcwd())
    logging.info(somalier_samples_fpath)
    df = pd.read_csv(somalier_samples_fpath, delimiter='\t')
    pairs_df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    with to_path(somalier_samples_fpath).open() as f:
        inferred_ped = Ped(f)
    with to_path(expected_ped_fpath).open() as f:
        expected_ped = Ped(f)

    bad = df.gt_depth_mean == 0.0
    if bad.any():
        warning(
            f'⚠️ Excluded {len(df[bad])}/{len(df)} samples with zero '
            f'mean GT depth from pedigree/sex checks: {", ".join(df[bad].sample_id)}',
        )
        info('')
    bad_ids = list(df[bad].sample_id)  # for checking in pairs_df
    df = df[~bad]

    info('*Inferred vs. reported sex:*')
    # Rename Ped sex to human-readable tags
    df.sex = df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    df.original_pedigree_sex = df.original_pedigree_sex.apply(lambda x: {'-9': 'unknown'}.get(x, x))
    missing_inferred_sex = df.sex == 'unknown'
    missing_provided_sex = df.original_pedigree_sex == 'unknown'
    mismatching_female = (df.sex == 'female') & (df.original_pedigree_sex == 'male')
    mismatching_male = (df.sex == 'male') & (df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male
    mismatching_other = (df.sex != df.original_pedigree_sex) & (~mismatching_female) & (~mismatching_male)
    matching_sex = ~mismatching_sex & ~mismatching_other

    def _print_stats(df_filter):
        for _, row_ in df[df_filter].iterrows():
            info(
                f' {row_.sample_id} ('
                f'provided: {row_.original_pedigree_sex}, '
                f'inferred: {row_.sex}, '
                f'mean depth: {row_.gt_depth_mean})',
            )

    if mismatching_sex.any():
        info(f'❗ {len(df[mismatching_sex])}/{len(df)} PED samples with mismatching sex:')
        _print_stats(mismatching_sex)
    if missing_provided_sex.any():
        info(f'⚠️ {len(df[missing_provided_sex])}/{len(df)} samples with missing provided sex:')
        _print_stats(missing_provided_sex)
    if missing_inferred_sex.any():
        info(f'⚠️ {len(df[missing_inferred_sex])}/{len(df)} samples with failed inferred sex:')
        _print_stats(missing_inferred_sex)
    inferred_cnt = len(df[~missing_inferred_sex])
    matching_cnt = len(df[matching_sex])
    info(
        f'✅ Sex inferred for {inferred_cnt}/{len(df)} samples, matching '
        f'for {matching_cnt if matching_cnt != inferred_cnt else "all"} samples.',
    )
    info('')

    info('*Relatedness:*')
    expected_ped_sample_by_id = {s.sample_id: s for s in expected_ped.samples()}
    inferred_ped_sample_by_id = {s.sample_id: s for s in inferred_ped.samples()}

    mismatching_unrelated_to_related = []
    mismatching_related_to_unrelated = []

    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        expected_ped_s1 = expected_ped_sample_by_id.get(s1)
        expected_ped_s2 = expected_ped_sample_by_id.get(s2)
        inferred_ped_s1 = inferred_ped_sample_by_id.get(s1)
        inferred_ped_s2 = inferred_ped_sample_by_id.get(s2)
        # Suppressing all logging output from peddy, otherwise it would clutter the logs
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            if expected_ped_s1 and expected_ped_s2:
                expected_rel = expected_ped.relation(expected_ped_s1, expected_ped_s2)
            else:
                expected_rel = 'unknown'
            if inferred_ped_s1 and inferred_ped_s2:
                inferred_rel = inferred_ped.relation(inferred_ped_s1, inferred_ped_s2)
            else:
                inferred_rel = 'unknown'

        if inferred_rel != expected_rel:
            # Constructing a line for a report:
            line = ''
            if (fam1 := expected_ped_s1.family_id if expected_ped_s1 else None) == (
                fam2 := expected_ped_s2.family_id if expected_ped_s2 else None
            ):
                line += f'{fam1}: {s1} - {s2}'
            else:
                line += s1 + (f' ({fam1})' if fam1 and fam1 != s1 else '')
                line += ' - '
                line += s2 + (f' ({fam2})' if fam2 and fam2 != s2 else '')
            line = (
                f'{line}, '
                f'provided: "{expected_rel}", '
                f'inferred: "{inferred_rel}", '
                f'kin={row["relatedness"]}, '
                f'ibs0={row["ibs0"]}, '
                f'ibs2={row["ibs2"]}'
            )

            if (
                expected_rel == 'unknown'
                and inferred_rel != 'unknown'
                or expected_rel == 'unrelated'
                and inferred_rel != 'unrelated'
            ):
                if row['relatedness'] > 0.1:
                    mismatching_unrelated_to_related.append(line)
            else:
                mismatching_related_to_unrelated.append(line)

        pairs_df.loc[idx, 'provided_rel'] = expected_rel  # type: ignore
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel  # type: ignore

    if mismatching_unrelated_to_related:
        info(
            f'⚠️ Found {len(mismatching_unrelated_to_related)} '
            f'sample pair(s) that are provided as unrelated, are inferred as '
            f'related:',
        )
        for i, pair in enumerate(mismatching_unrelated_to_related):
            info(f' {i + 1}. {pair}')
    if mismatching_related_to_unrelated:
        info(
            f'❗ Found {len(mismatching_related_to_unrelated)} sample pair(s) '
            f'that are provided as related, but inferred as unrelated:',
        )
        for i, pair in enumerate(mismatching_related_to_unrelated):
            info(f' {i + 1}. {pair}')
    if not mismatching_unrelated_to_related and not mismatching_related_to_unrelated:
        info('✅ Inferred pedigree matches for all provided related pairs.')
    info('')

    print_contents(
        df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
    )

    # Constructing Slack message
    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "Somalier pedigree report"}>'
    elif not title:
        title = 'Somalier pedigree report'
    text = '\n'.join([title] + _messages)

    if send_to_slack:
        send_message(text)


def print_contents(
    samples_df,
    pairs_df,
    somalier_samples_fpath,
    somalier_pairs_fpath,
):
    """
    Print useful information to manually review pedigree check results
    """
    if len(samples_df) < 400:
        samples_str = samples_df.to_string()
        logging.info(f'Somalier results, samples (based on {somalier_samples_fpath}):\n{samples_str}\n')
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
        logging.info(f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n{pairs_str}\n')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
