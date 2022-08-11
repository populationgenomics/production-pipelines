#!/usr/bin/env python3

"""
Checks metrics in MultiQC output, based on thresholds in the `qc_thresholds` 
config section.

Script can send a report to a Slack channel. To enable that, set `SLACK_TOKEN`
and `SLACK_CHANNEL` environment variables, and add "Seqr Loader" app into 
a channel with:

/invite @Seqr Loader
"""

import logging
import json
import os
from typing import Optional

import click
from cpg_utils import to_path
from cpg_utils.config import get_config

logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


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
    sequencing_type = get_config()['workflow']['sequencing_type']
    thresholds_d = get_config()['qc_thresholds'].get(sequencing_type)
    min_thresholds_d = thresholds_d.get('min')
    max_thresholds_d = thresholds_d.get('max')

    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
    lines = []
    for tool_d in d['report_general_stats_data']:
        for sample, failed_val_by_metric in tool_d.items():
            for metric, min_value in min_thresholds_d.items():
                if metric in failed_val_by_metric:
                    val = failed_val_by_metric[metric]
                    if val < min_value:
                        val = float(val)
                        if isinstance(val, float):
                            val_str = f'{0:.2g}'.format(val)
                        else:
                            val_str = str(val)
                        lines.append(f'{sample}: {metric}={val_str} (< {min_value})')

            for metric, min_value in max_thresholds_d.items():
                if metric in failed_val_by_metric:
                    val = failed_val_by_metric[metric]
                    if val > min_value:
                        if isinstance(val, float):
                            val_str = f'{0:.2g}'.format(val)
                        else:
                            val_str = str(val)
                        lines.append(f'{sample}: {metric}={val_str} (> {min_value})')

    slack_channel = os.environ.get('SLACK_CHANNEL')
    slack_token = os.environ.get('SLACK_TOKEN')
    if send_to_slack and slack_token and slack_channel:
        from slack_sdk.errors import SlackApiError
        from slack_sdk import WebClient

        slack_client = WebClient(token=slack_token)
        try:
            title = title or 'MutliQC report'
            text = f'*[{dataset}]* <{html_url}|{title}>\n'
            text += '\n'.join(lines)

            slack_client.api_call(  # pylint: disable=duplicate-code
                'chat.postMessage',
                json={
                    'channel': slack_channel,
                    'text': text,
                },
            )
        except SlackApiError as err:
            logging.error(f'Error posting to Slack: {err}')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
