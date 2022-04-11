"""
Abstract reporter of stages statuses.
"""

import logging
from abc import ABC, abstractmethod
from enum import Enum

from hailtop.batch.job import Job
from hailtop.batch import Batch

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
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums
    /analysis.py#L14-L21
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

    @abstractmethod
    def add_updaters_jobs(
        self,
        b: Batch,
        output: Path | dict[str, Path],
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
    ):
        """
        Add Hail Batch jobs that update the analysis status.
        """
