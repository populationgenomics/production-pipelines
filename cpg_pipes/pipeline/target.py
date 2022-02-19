"""
Defines a target that stage can act upon. Classes like Sample, Project, Pipeline
expend this class.
"""
from abc import abstractmethod
from typing import Dict

from cpg_pipes.smdb.types import AnalysisType, Analysis


class Target:
    """
    Stage target: e.g. a sample, a dataset/project, entire cohort/pipeline.
    """
    def __init__(self):
        # From SMDB Analysis entries:
        self.analysis_by_type: Dict[AnalysisType, Analysis] = dict()
        # Whether to process even if outputs exist:
        self.forced = False
        # If not set, exclude from the pipeline:
        self.active = True

    @property
    @abstractmethod
    def unique_id(self) -> str:
        """
        ID should be unique across target of all levels.
        """
        pass
