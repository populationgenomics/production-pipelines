"""
Defines a target that stage can act upon. Classes like Sample, Project, Pipeline
expend this class.
"""
from dataclasses import dataclass, field, KW_ONLY
from typing import Dict
from cpg_pipes.smdb.types import AnalysisType, Analysis


@dataclass
class Target:
    """
    Stage target: e.g. a sample, a dataset/project, entire cohort/pipeline.
    """

    _: KW_ONLY
    # From SMDB Analysis entries:
    analysis_by_type: Dict[AnalysisType, Analysis] = field(
        default_factory=dict, repr=False
    )
    # Whether to process even if outputs exist:
    forced = False
    # If not set, exclude from the pipeline:
    active = True

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
