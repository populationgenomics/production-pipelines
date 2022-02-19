import functools
import inspect
import logging
import sys
from typing import List, Optional, Callable, Type, Union

from cpg_pipes.pipeline import Pipeline
from cpg_pipes.pipeline.stage import Stage
from cpg_pipes.smdb.types import AnalysisType

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


StageDecorator = Callable[..., 'Stage']


def stage(
    _cls: Optional[Type[Stage]] = None, 
    *,
    sm_analysis_type: Optional[AnalysisType] = None, 
    requires_stages: Optional[Union[List[StageDecorator], StageDecorator]] = None,
    skipped: bool = False,
    required: bool = True,
    assume_results_exist: bool = False,
    forced: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with an optional argument.
    The goal is to allow cleaner defining of custom pipeline stages, without
    requiring to implement constructor. E.g.

    @stage(sm_analysis_type=AnalysisType.GVCF, requires_stages=CramStage)
    class GvcfStage(SampleStage):
        def expected_result(self, sample: Sample):
            ...
        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
            ...
    """
    def decorator_stage(cls) -> StageDecorator:
        @functools.wraps(cls)
        def wrapper_stage(pipeline: 'Pipeline') -> Stage:
            return cls(
                name=cls.__name__,
                pipeline=pipeline,
                requires_stages=requires_stages,
                sm_analysis_type=sm_analysis_type,
                skipped=skipped,
                required=required,
                assume_results_exist=assume_results_exist, 
                forced=forced,
            )
        return wrapper_stage

    if _cls is None:
        return decorator_stage
    else:
        return decorator_stage(_cls)


def skipped(
    _fun: Optional[StageDecorator] = None, 
    *,
    assume_results_exist: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Decorator on top of `@stage` that sets the self.skipped field to True

    @skipped
    @stage
    class MyStage1(SampleStage):
        ...

    @skipped
    @stage(assume_results_exist=True)
    class MyStage2(SampleStage):
        ...
    """
    def decorator_stage(fun) -> StageDecorator:
        @functools.wraps(fun)
        def wrapper_stage(*args, **kwargs) -> Stage:
            stage = fun(*args, **kwargs)
            stage.skipped = True
            stage.assume_results_exist = assume_results_exist
            return stage

        return wrapper_stage

    if _fun is None:
        return decorator_stage
    else:
        return decorator_stage(_fun)


def find_stages_in_module(module_name: str) -> List[StageDecorator]:
    """
    Get all declared stages in a module by its name
    >>> pipeline = Pipeline()
    >>> pipeline.set_stages(find_stages_in_module(__name__))
    """
    stages = []
    for obj in vars(sys.modules[module_name]).values():
        # @stage decorators are functions
        if not inspect.isfunction(obj): 
            continue
        # decorators add "return" annotations
        if 'return' not in obj.__annotations__: 
            continue
        try:
            ret_name = obj.__annotations__['return'].__name__
        except AttributeError:
            continue
        if ret_name == 'Stage':
            stages.append(obj)
    
    return stages
