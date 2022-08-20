#!/usr/bin/env python3

"""
Checks metrics in MultiQC output, based on thresholds in the `qc_thresholds` 
config section.

Script can send a report to a Slack channel. To enable that, set `SLACK_TOKEN`
and `SLACK_CHANNEL` environment variables, and add "Seqr Loader" app into 
a channel with:

/invite @Seqr Loader
"""
import dataclasses
import logging
import json
import pprint
import os
from collections import defaultdict
from typing import Optional, Literal

import click
from cpg_utils import to_path
from cpg_utils.config import get_config

logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(get_config()['workflow'].get('log_level', 'INFO'))


@click.command()
@click.option(
    '--multiqc-json',
    'multiqc_json_path',
    required=True,
    help='Path to MultiQC JSON output',
)
@click.option(
    '--html-url',
    'html_url',
    required=True,
    help='MultiQC HTML URL',
)
@click.option('--dataset', 'dataset', help='Dataset name')
@click.option('--title', 'title', help='Report title')
@click.option(
    '--send-to-slack/--no-send-to-slack',
    'send_to_slack',
    help='Send log to Slack message, according to environment variables SLACK_CHANNEL and SLACK_TOKEN',
)
def main(
    multiqc_json_path: str,
    html_url: str,
    dataset: str,
    title: Optional[str] = None,
    send_to_slack: bool = True,
):
    """
    Check metrics in MultiQC json and send info about failed samples
    as a Slack message.
    """
    seq_type = get_config()['workflow']['sequencing_type']

    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
        sections = d['report_general_stats_data']
        print(f'report_general_stats_data: {pprint.pformat(sections)}')

    bad_lines_by_sample = defaultdict(list)
    for config_key, fail_sign, good_sign, is_fail in [
        (
            'min',
            '<',
            '≥',
            lambda val_, thresh_: val_ < thresh_,
        ),
        (
            'max',
            '>',
            '≤',
            lambda val_, thresh_: val_ > thresh_,
        ),
    ]:
        threshold_d = (
            get_config()['qc_thresholds'].get(seq_type, {}).get(config_key, {})
        )
        for section in sections:
            for sample, val_by_metric in section.items():
                for metric, threshold in threshold_d.items():
                    if metric in val_by_metric:
                        val = val_by_metric[metric]
                        if is_fail(val, threshold):
                            line = f'{metric}={val:0.2f}{fail_sign}{threshold:0.2f}'
                            bad_lines_by_sample[sample].append(line)
                            logger.debug(f'⭕ {sample}: {line}')
                        else:
                            line = f'{metric}={val:0.2f}{good_sign}{threshold:0.2f}'
                            logger.debug(f'✅ {sample}: {line}')

    # Constructing Slack message
    title = f'*[{dataset}]* <{html_url}|{title or "MutliQC report"}>'
    lines = []
    if bad_lines_by_sample:
        lines.append(f'{title}')
        for sample, bad_lines in bad_lines_by_sample.items():
            lines.append(f'⭕ {sample}: ' + ', '.join(bad_lines))
    else:
        lines.append(f'✅ {title}')
    message = '\n'.join(lines)

    slack_channel = get_config().get('slack', {}).get('channel')
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
                    'text': message,
                },
            )
        except SlackApiError as err:
            logging.error(f'Error posting to Slack: {err}')
    else:
        print(message)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
