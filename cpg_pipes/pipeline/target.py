"""
Defines a target that stage can act upon. Classes like Sample, Dataset, Pipeline
extend this class.
"""

from typing import Dict
from cpg_pipes.smdb.types import AnalysisType, Analysis


class Target:
    """
    Stage for a target
    """
    def __init__(self):
        # From SMDB Analysis entries:
        self.analysis_by_type: Dict[AnalysisType, Analysis] = dict()
        # Whether to process even if outputs exist:
        self.forced: bool = False
        # If not set, exclude from the pipeline:
        self.active: bool = True

    @property
    def unique_id(self) -> str:
        """
        ID should be unique across target of all levels. We are raising
        NotImplementedError instead of making it abstractclass because
        mypy is not happy about binding TypeVar to abstract classes:
        https://stackoverflow.com/questions/48349054/how-do-you-annotate-the-type-of-an-abstract-class-with-mypy
        Specifically,
        ```
        TypeVar('TargetT', bound=Target)
        ```
        Will raise:
        ```
        Only concrete class can be given where "Type[Target]" is expected
        ```
        
        """
        raise NotImplementedError
