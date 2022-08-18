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

    @abstractmethod
    def add_updaters_jobs(
        self,
        b: Batch,
        output: str | Path | Resource | dict[str, Path | Resource],
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
    ) -> list[Job]:
        """
        Add Hail Batch jobs that update the analysis status.
        """

    @abstractmethod
    def create_analysis(
        self,
        output: str,
        analysis_type: str,
        analysis_status: str,
        target: Target,
        meta: dict | None = None,
        project_name: str = None,
    ) -> int | None:
        """Record analysis entry"""
