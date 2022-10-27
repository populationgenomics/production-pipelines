"""
Sending notifications to Slack. To enable, create a channel, add "Seqr Loader" app
into a channel with:

/invite @Seqr Loader

Make sure `slack/channel`, `slack/token_secret_id`, and `slack/token_project_id`
configuration values are set.
"""

from textwrap import dedent

from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config

from hailtop.batch.job import Job


def get_token() -> str:
    """
    Returns Slack token.
    """
    token_secret_id = get_config()['slack'].get('token_secret_id')
    token_project_id = get_config()['slack'].get('token_project_id')
    if not token_secret_id or not token_project_id:
        raise ValueError(
            '`slack.token_secret_id` and `slack.token_project_id` '
            'must be set in config to retrieve Slack token'
        )
    token = read_secret(
        project_id=token_project_id,
        secret_name=token_secret_id,
        fail_gracefully=False,
    )
    assert token
    return token


def slack_env(j: Job):
    """
    Add environment variables that configure Slack reporter.
    """
    if not (channel := get_config().get('slack', {}).get('channel')):
        return None
    j.env('SLACK_CHANNEL', channel)
    j.env('SLACK_TOKEN', get_token())


def slack_message_cmd(
    j: Job,
    text: str | None = None,
    data: dict[str, str] | None = None,
):
    """
    Make a bash command that prepares and sends a Slack message.
    Message can be constructed from a dict `data`, or can be passed directly
    as text (`text`). In either case, strings can use Slack `mrkdwn`
    for formatting: https://api.slack.com/reference/surfaces/formatting

    Assumes `slack/channel` configuration parameter is set, as well as
    SLACK_TOKEN environment variable.
    """
    if not (channel := get_config()['slack'].get('channel')):
        return ''
    msg = ''
    if data:
        msg += '\\n'.join(f'{k}: {v}' for k, v in data.items())
    if text:
        msg += f'\\n{text}'
    if not msg:
        return ''
    slack_env(j)
    j.command(
        dedent(
            f"""
    curl -X POST \
    -H "Authorization: Bearer {get_token()}" \
    -H "Content-type: application/json" \
    -d '{{"channel": "'{channel}'", "text": "{msg}"}}' \
    https://slack.com/api/chat.postMessage
    """
        )
    )
