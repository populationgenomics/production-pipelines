#!/usr/bin/env python3

"""
This script parses "somalier relate" (https://github.com/brentp/somalier) outputs,
and returns a report whether sex and pedigree matches the provided PED file.

You can have report sent to channel. To do that, add "Seqr Loader" app into a channel:

/invite @Seqr Loader

And run the script with "--slack-channel channel_name".
"""

import contextlib
import logging
import os
from io import StringIO
from typing import Optional, Tuple, List

import click
import pandas as pd
import slack_sdk
from peddy import Ped
from slack_sdk.errors import SlackApiError

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
@click.option('--dataset', 'dataset', help='Dataset name')
def main(
    somalier_samples_fpath: str,
    somalier_pairs_fpath: str,
    dataset: Optional[str] = None,
):
    """
    Report pedigree inconsistencies, given somalier outputs.
    """
    pedlog = PedigreeLogger()
    check_pedigree(
        somalier_samples_fpath,
        somalier_pairs_fpath,
        dataset=dataset,
        pedlog=pedlog,
    )

    slack_token = os.environ.get('SLACK_TOKEN')
    slack_channel = os.environ.get('SLACK_CHANNEL')
    if slack_token and slack_channel:
        slack_client = slack_sdk.WebClient(token=slack_token)
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
    fp = StringIO()
    df.to_csv(fp, sep='\t', index=False)
    ped = Ped(StringIO(fp.getvalue()))

    pedlog = pedlog or PedigreeLogger()

    bad = df.gt_depth_mean == 0.0
    if bad.any():
        pedlog.warning(
            f'*[{dataset}]* excluding {len(df[bad])}/{len(df)} samples with zero '
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
                f'\t{row_.sample_id} ('
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
        f'in {matching_cnt if matching_cnt != inferred_cnt else "all"} samples.'
    )
    pedlog.info('')

    pedlog.info(f'*[{dataset}] relatedness:*')
    sample_by_id = {s.sample_id: s for s in ped.samples()}

    pairs_provided_as_unrelated_but_inferred_related = []
    other_mismatching_pairs = []

    for idx, row in pairs_df.iterrows():
        s1 = row['#sample_a']
        s2 = row['sample_b']
        if s1 in bad_ids or s2 in bad_ids:
            continue

        inferred_rel, reason = infer_relationship(
            row['relatedness'], row['ibs0'], row['ibs2']
        )
        # Supressing all logging output from peddy, otherwise it would clutter the logs
        with contextlib.redirect_stderr(None), contextlib.redirect_stdout(None):
            peddy_rel = ped.relation(sample_by_id[s1], sample_by_id[s2])

        def _match_peddy_with_inferred(peddy_rel_):
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
            }.get(peddy_rel_)

        if (
            _match_peddy_with_inferred(peddy_rel) == 'unknown'
            and inferred_rel != 'unknown'
            or _match_peddy_with_inferred(peddy_rel) == 'unrelated'
            and inferred_rel != 'unrelated'
        ):
            pairs_provided_as_unrelated_but_inferred_related.append(
                f'"{s1}" and "{s2}", '
                f'provided: "{peddy_rel}", '
                f'inferred: "{inferred_rel}" ({reason})'
            )
        elif inferred_rel != _match_peddy_with_inferred(peddy_rel):
            other_mismatching_pairs.append(
                f'"{s1}" and "{s2}", '
                f'provided: "{peddy_rel}", '
                f'inferred: "{inferred_rel}" ({reason})'
            )
        pairs_df.loc[idx, 'provided_rel'] = peddy_rel
        pairs_df.loc[idx, 'inferred_rel'] = inferred_rel

    if pairs_provided_as_unrelated_but_inferred_related:
        pedlog.info(
            f'{len(pairs_provided_as_unrelated_but_inferred_related)} sample pair(s) '
            f'that are provided as unrelated, but inferred as related:'
        )
        for i, pair in enumerate(pairs_provided_as_unrelated_but_inferred_related):
            pedlog.info(f'#{i + 1} {pair}')
    if other_mismatching_pairs:
        pedlog.info(
            f'{len(other_mismatching_pairs)} sample pair(s) with mismatched relatedness:'
        )
        for i, pair in enumerate(other_mismatching_pairs):
            pedlog.info(f'#{i + 1} {pair}')
    else:
        pedlog.info(f'Inferred pedigree matches for all {len(df)} samples.')
    pedlog.info('')

    print_contents(
        df,
        pairs_df,
        somalier_samples_fpath,
        somalier_pairs_fpath,
    )

    pedlog.mismatching_sex = mismatching_sex.any()
    pedlog.mismatching_relationships = other_mismatching_pairs
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
                reason + f' ibs0/ibs2={(ibs0 / ibs2):.4f} < 0.005',
            )
        elif 0.015 < (ibs0 / ibs2) < 0.052:
            result = (
                'siblings',
                reason + f' 0.015 < ibs0/ibs2={(ibs0 / ibs2):.4f} < 0.052',
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
