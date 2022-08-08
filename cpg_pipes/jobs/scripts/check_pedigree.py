#!/usr/bin/env python3

"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

Script can send a report to a Slack channel. To enable that, set `SLACK_TOKEN`
and `SLACK_CHANNEL` environment variables, and add "Seqr Loader" app into 
a channel with:

/invite @Seqr Loader
"""

import contextlib
import logging
import os
from typing import Optional, Tuple, List

import click
import pandas as pd
from cpg_utils import to_path
from peddy import Ped

logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger = logging.getLogger(__file__)
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
    '--expected-ped',
    'expected_ped_fpath',
    required=True,
    help='Path to PED file with expected pedigree',
)
@click.option('--dataset', 'dataset', help='Dataset name')
@click.option(
    '--send-to-slack/--no-send-to-slack',
    'send_to_slack',
    help='Send log to Slack message, according to environment variables SLACK_CHANNEL and SLACK_TOKEN',
)
def main(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    dataset: Optional[str] = None,
    send_to_slack: bool = True,
):
    """
    Report pedigree inconsistencies, given somalier outputs.
    """
    pedlog = PedigreeLogger()
    check_pedigree(
        somalier_samples_fpath,
        somalier_pairs_fpath,
        expected_ped_fpath,
        dataset=dataset,
        pedlog=pedlog,
    )

    slack_channel = os.environ.get('SLACK_CHANNEL')
    slack_token = os.environ.get('SLACK_TOKEN')
    if send_to_slack and slack_token and slack_channel:
        from slack_sdk.errors import SlackApiError
        from slack_sdk import WebClient

        slack_client = WebClient(token=slack_token)
        try:
            slack_client.api_call(  # pylint: disable=duplicate-code
                'chat.postMessage',
                json={
                    'channel': slack_channel,
                    'text': '\n'.join(pedlog.messages),
                },
            )
        except SlackApiError as err:
            logging.error(f'Error posting to Slack: {err}')


class PedigreeLogger:
    """
    Recording messages to send a report to Slack.
    """

    def __init__(self):
        self.mismatching_sex = False
        self.mismatching_relationships = False
        self.messages: List[str] = []

    def info(self, msg):
        """
        Record and forward.
        """
        self.messages.append(msg)
        logger.info(msg)

    def warning(self, msg):
        """
        Record and forward.
        """
        self.messages.append(msg)
        logger.info(msg)

    def error(self, msg):
        """
        Record and forward.
        """
        self.messages.append(msg)
        logger.error(msg)


def check_pedigree(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    expected_ped_fpath: str,
    dataset: Optional[str] = None,
    pedlog: Optional[PedigreeLogger] = None,
) -> PedigreeLogger:
    """
    Report pedigree inconsistencies, given somalier outputs.
    """
    print(os.getcwd())
    print(somalier_samples_fpath)
    df = pd.read_csv(somalier_samples_fpath, delimiter='\t')
    pairs_df = pd.read_csv(somalier_pairs_fpath, delimiter='\t')
    with to_path(expected_ped_fpath).open() as f:
        expected_ped = Ped(f)

    pedlog = pedlog or PedigreeLogger()

    bad = df.gt_depth_mean == 0.0
    if bad.any():
        pedlog.warning(
            f'*[{dataset}]* excluded {len(df[bad])}/{len(df)} samples with zero '
            f'mean GT depth from pedigree/sex checks: {", ".join(df[bad].sample_id)}'
        )
        pedlog.info('')
    bad_ids = list(df[bad].sample_id)  # for checking in pairs_df
    df = df[~bad]

    pedlog.info(f'*[{dataset}] inferred vs. reported sex:*')
    # Rename Ped sex to human-readable tags
    df.sex = df.sex.apply(lambda x: {1: 'male', 2: 'female'}.get(x, 'unknown'))
    df.original_pedigree_sex = df.original_pedigree_sex.apply(
        lambda x: {'-9': 'unknown'}.get(x, x)
    )
    missing_inferred_sex = df.sex == 'unknown'
    missing_provided_sex = df.original_pedigree_sex == 'unknown'
    mismatching_female = (df.sex == 'female') & (df.original_pedigree_sex == 'male')
    mismatching_male = (df.sex == 'male') & (df.original_pedigree_sex == 'female')
    mismatching_sex = mismatching_female | mismatching_male
    mismatching_other = (
        (df.sex != df.original_pedigree_sex)
        & (~mismatching_female)
        & (~mismatching_male)
    )
    matching_sex = ~mismatching_sex & ~mismatching_other

    def _print_stats(df_filter):
        for _, row_ in df[df_filter].iterrows():
            pedlog.info(
                f' {row_.sample_id} ('
                f'provided: {row_.original_pedigree_sex}, '
                f'inferred: {row_.sex}, '
                f'mean depth: {row_.gt_depth_mean})'
            )

    if mismatching_sex.any():
        pedlog.info(
            f'{len(df[mismatching_sex])}/{len(df)} PED samples with mismatching sex:'
        )
        _print_stats(mismatching_sex)
    if missing_provided_sex.any():
        pedlog.info(
            f'{len(df[missing_provided_sex])}/{len(df)} samples with missing provided sex:'
        )
        _print_stats(missing_provided_sex)
    if missing_inferred_sex.any():
        pedlog.info(
            f'{len(df[missing_inferred_sex])}/{len(df)} samples with failed inferred sex:'
        )
        _print_stats(missing_inferred_sex)
    inferred_cnt = len(df[~missing_inferred_sex])
    matching_cnt = len(df[matching_sex])
    pedlog.info(
        f'Sex inferred for {inferred_cnt}/{len(df)} samples, matching '
        f'for {matching_cnt if matching_cnt != inferred_cnt else "all"} samples.'
    )
    pedlog.info('')

    pedlog.info(f'*[{dataset}] relatedness:*')
    ped_sample_by_id = {s.sample_id: s for s in expected_ped.samples()}

    mismatching_unrelated_to_related = []
    mismatching_unrelated_to_closely_related = []
    mismatching_related_to_unrelated = []

    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        inferred_rel, reason = infer_relationship(
            row['relatedness'], row['ibs0'], row['ibs2']
        )

        ped_s1 = ped_sample_by_id.get('s1')
        ped_s2 = ped_sample_by_id.get('s2')
        if ped_s1 and ped_s2:
            # Supressing all logging output from peddy, otherwise it would clutter the logs
            with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
                peddy_rel = expected_ped.relation(ped_s1, ped_s2)
        else:
            peddy_rel = 'unknown'
        matched_peddy_rel = {
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

        def _repr_cur_pair() -> str:
            fam1 = expected_ped.get(sample_id=s1).family_id
            fam2 = expected_ped.get(sample_id=s2).family_id
            line = ''
            if fam1 == fam2:
                line += f'{fam1}: {s1} - {s2}'
            else:
                line += s2 + (f' ({fam1})' if fam1 != s1 else '')
                line += ' - '
                line += s2 + (f' ({fam2})' if fam2 != s2 else '')
            return line + (
                f', provided: "{peddy_rel}", inferred: "{inferred_rel}" ({reason})'
            )

        if (
            matched_peddy_rel == 'unknown'
            and inferred_rel != 'unknown'
            or matched_peddy_rel == 'unrelated'
            and inferred_rel != 'unrelated'
        ):
            mismatching_unrelated_to_related.append(_repr_cur_pair())
            if is_close(row['relatedness']):
                mismatching_unrelated_to_closely_related.append(_repr_cur_pair())

        elif inferred_rel != matched_peddy_rel:
            mismatching_related_to_unrelated.append(_repr_cur_pair())

        pairs_df.loc[idx, 'provided_rel'] = peddy_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    if mismatching_unrelated_to_related:
        pedlog.info(
            f'Found {len(mismatching_unrelated_to_related)} sample '
            f'pair(s) that are provided as unrelated, but inferred as related '
            f'below first degree'
        )

    if mismatching_unrelated_to_closely_related:
        pedlog.info(
            f'Found {len(mismatching_unrelated_to_closely_related)} '
            f'sample pair(s) that are provided as unrelated, are inferred as '
            f'first-degree related, twins or identical:'
        )
        for i, pair in enumerate(mismatching_unrelated_to_closely_related):
            pedlog.info(f' {i + 1}. {pair}')
    if mismatching_related_to_unrelated:
        pedlog.info(
            f'Found {len(mismatching_related_to_unrelated)} sample pair(s) '
            f'that are provided as related, but inferred as unrelated:'
        )
        for i, pair in enumerate(mismatching_related_to_unrelated):
            pedlog.info(f' {i + 1}. {pair}')
    if (
        not mismatching_unrelated_to_related
        and not mismatching_related_to_unrelated
        and not mismatching_related_to_unrelated
    ):
        pedlog.info(f'Inferred pedigree matches for all provided related pairs.')
    pedlog.info('')

    print_contents(
        df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
    )

    pedlog.mismatching_sex = mismatching_sex.any()
    pedlog.mismatching_relationships = mismatching_related_to_unrelated
    return pedlog


def infer_relationship(kin: float, ibs0: float, ibs2: float) -> Tuple[str, str]:
    """
    Inferres relashionship labels based on the kin coefficient
    and ibs0 and ibs2 values. Returns relashionship label and reason.
    """
    if kin < 0.1:
        result = ('unrelated', f'kin={kin:.3f} < 0.1')
    elif kin < 0.38:
        result = ('below_first_degree', f'kin={kin:.3f} < 0.38')
    elif kin <= 0.62:
        reason = f'kin={kin:.3f} < 0.62'
        if (ibs0 / ibs2) < 0.005:
            result = (
                'parent-child',
                reason + f', ibs0/ibs2={(ibs0 / ibs2):.4f} < 0.005',
            )
        elif 0.015 < (ibs0 / ibs2) < 0.052:
            result = (
                'siblings',
                reason + f', 0.015 < ibs0/ibs2={(ibs0 / ibs2):.4f} < 0.052',
            )
        else:
            result = (
                'first_degree',
                reason
                + (
                    f', not parent-child (ibs0/ibs2={(ibs0 / ibs2):.4f} > 0.005)'
                    f', not sibling (0.015 < ibs0/ibs2={(ibs0 / ibs2):.4f} < 0.052)'
                ),
            )
    elif kin < 0.8:
        result = ('first_degree_or_duplicate_or_twins', f'kin={kin:.3f} < 0.8')
    else:
        assert kin >= 0.8
        result = ('duplicate_or_twins', f'kin={kin:.3f} >= 0.8')
    return result


def is_close(kin: float) -> bool:
    """First degree, duplicate or twin"""
    return kin > 0.62


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
        logger.info(
            f'Somalier results, samples (based on {somalier_samples_fpath}):\n'
            f'{samples_str}\n'
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
            f'Somalier results, sample pairs (based on {somalier_pairs_fpath}):\n'
            f'{pairs_str}\n'
        )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
