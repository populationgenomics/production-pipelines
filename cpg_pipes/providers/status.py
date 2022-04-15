"""
Abstract reporter of stages statuses.
"""

import logging
from abc import ABC, abstractmethod
from enum import Enum

from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from .. import Path
from ..targets import Target

logger = logging.getLogger(__file__)


class StatusReporterError(Exception):
    """
    Raised for problems reporting status
    """


class AnalysisStatus(Enum):
    """
    Corresponds to SMDB Analysis statuses:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L14-L21
    """

    QUEUED = 'queued'
    IN_PROGRESS = 'in-progress'
    FAILED = 'failed'
    COMPLETED = 'completed'
    UNKNOWN = 'unknown'

    @staticmethod
    def parse(name: str) -> 'AnalysisStatus':
        """
        Parse str and create a AnalysisStatus object
        """
        return {v.value: v for v in AnalysisStatus}[name.lower()]


class StatusReporter(ABC):
    """
    Status reporter
    """
    def __init__(self):
        self.slack_channel = None
        self.slack_token = None

    @abstractmethod
    def add_updaters_jobs(
        self,
        b: Batch,
        output: Path | Resource | dict[str, Path | Resource],
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
    ):
        """
        Add Hail Batch jobs that update the analysis status.
        """

    def slack_env(self, j: Job):
        """
        Add environment variables that configure Slack reporter.
        """
        if self.slack_channel and self.slack_token:
            j.env('SLACK_CHANNEL', self.slack_channel)
            j.env('SLACK_TOKEN', self.slack_token)

    def slack_message_cmd(
        self, 
        text: str | None = None, 
        data: dict[str, str] | None = None,
    ) -> str:
        """
        Add command to the job that sends Slack message.
        Message can be a dictionary (`data`), which will be correspondingly formatted;
        or text (`text`). Either can use Slack mrkdwn for formatting: 
        https://api.slack.com/reference/surfaces/formatting
        """
        if not self.slack_channel or not self.slack_token:
            return ''
        msg = ''
        if data:
            msg += '\\n'.join(f'{k}: {v}' for k, v in data.items())
        if text:
            msg += f'\\n{text}'
        if not msg:
            return ''
        return f"""
        curl -X POST \
        -H "Authorization: Bearer $SLACK_TOKEN" \
        -H "Content-type: application/json" \
        -d '{{"channel": "'$SLACK_CHANNEL'", "text": "{msg}"}}' \
        https://slack.com/api/chat.postMessage
        """
