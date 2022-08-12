"""
Provides a `Pipeline` class and a `@stage` decorator that allow to define 
pipeline stages and plug them together.

A `Stage` describes the job that are added to Hail Batch, and outputs that are expected
to be produced. Each stage acts on a `Target`, which can be of a different level:
a `Sample`, a `Dataset`, a `Cohort` (= all input datasets combined). Pipeline plugs
stages together by resolving dependencies between different levels accordingly.

Examples of pipelines can be found in the `pipelines/` folder in the repository root.
"""
import functools
import logging
import pathlib
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from typing import cast, Callable, Union, TypeVar, Generic, Any, Optional, Type

from cloudpathlib import CloudPath
import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils.config import get_config
from cpg_utils import Path, to_path

from cpg_pipes.hb.batch import setup_batch, RegisteringBatch
from cpg_pipes.status import MetamistStatusReporter
from cpg_pipes.targets import Target, Dataset, Sample, Cohort
from cpg_pipes.status import StatusReporter
from cpg_pipes.utils import exists, timestamp, slugify
from cpg_pipes.inputs import get_cohort
from .exceptions import PipelineError, StageInputNotFound

logger = logging.getLogger(__file__)


# We record all initialised Stage subclasses, which we then use as a default
# list of stages when the user didn't pass them explicitly.
_ALL_DECLARED_STAGES = []


StageDecorator = Callable[..., 'Stage']

# Type variable to use with Generic to make sure a Stage subclass always matches the
# corresponding Target subclass. We can't just use the Target superclass because
# it would violate the Liskov substitution principle (i.e. any Stage subclass would
# have to be able to work on any Target subclass).
TargetT = TypeVar('TargetT', bound=Target)

ExpectedResultT = Union[Path, dict[str, Path], None]

StageOutputData = Union[Path, hb.Resource, dict[str, Path], dict[str, hb.Resource]]


# noinspection PyShadowingNames
class StageOutput:
    """
    Represents a result of a specific stage, which was run on a specific target.
    Can be a file path, or a Hail Batch Resource. Optionally wrapped in a dict.
    """

    def __init__(
        self,
        target: 'Target',
        data: StageOutputData | str | dict[str, str] | None = None,
        jobs: list[Job] | Job | None = None,
        meta: dict | None = None,
        reusable: bool = False,
        skipped: bool = False,
        error_msg: str | None = None,
        stage: Optional['Stage'] = None,
    ):
        # Converting str into Path objects.
        self.data: StageOutputData | None
        if isinstance(data, dict):
            self.data = {k: to_path(v) for k, v in data.items()}
        elif data is not None:
            self.data = to_path(data)
        else:
            self.data = data

        self.stage = stage
        self.target = target
        self.jobs: list[Job] = [jobs] if isinstance(jobs, Job) else (jobs or [])
        self.meta: dict = meta or {}
        self.reusable = reusable
        self.skipped = skipped
        self.error_msg = error_msg

    def __repr__(self) -> str:
        res = (
            f'StageOutput({self.data}'
            f' target={self.target}'
            f' stage={self.stage}'
            + (f' [reusable]' if self.reusable else '')
            + (f' [skipped]' if self.skipped else '')
            + (f' [error: {self.error_msg}]' if self.error_msg else '')
            + f' meta={self.meta}'
            + f')'
        )
        return res

    def as_path_or_resource(self, id=None) -> Path | hb.Resource:
        """
        Cast the result to Union[str, hb.Resource], error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        if self.data is None:
            raise ValueError(f'{self.stage}: output data is not available')

        if id is not None:
            if not isinstance(self.data, dict):
                raise ValueError(
                    f'{self.stage}: {self.data} is not a dictionary, can\'t get "{id}"'
                )
            return cast(dict, self.data)[id]

        if isinstance(self.data, dict):
            res = cast(dict, self.data)
            if len(res.values()) > 1:
                raise ValueError(
                    f'{res} is a dictionary with more than 1 element, '
                    f'please set the `id` parameter'
                )
            return list(res.values())[0]

        return self.data

    def as_path(self, id=None) -> Path:
        """
        Cast the result to path. Though exception if failed to cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, CloudPath | pathlib.Path):
            raise ValueError(f'{res} is not a path.')
        return cast(Path, res)

    def as_resource(self, id=None) -> hb.Resource:
        """
        Cast the result to Hail Batch Resource, or throw an error if the cast failed.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, hb.Resource):
            raise ValueError(f'{res} is not a Hail Batch Resource.')
        return cast(hb.Resource, res)

    def as_dict(self) -> dict[str, Path | hb.Resource]:
        """
        Cast the result to a dictionary, or throw an error if the cast failed.
        """
        if not isinstance(self.data, dict):
            raise ValueError(f'{self.data} is not a dictionary.')
        return self.data

    def as_resource_dict(self) -> dict[str, hb.Resource]:
        """
        Cast the result to a dictionary of Hail Batch Resources,
        or throw an error if the cast failed
        """
        return {k: self.as_resource(id=k) for k in self.as_dict()}

    def as_path_dict(self) -> dict[str, Path]:
        """
        Cast the result to a dictionary of strings,
        or throw an error if the cast failed.
        """
        return {k: self.as_path(id=k) for k in self.as_dict()}


# noinspection PyShadowingNames
class StageInput:
    """
    Represents an input for a stage run. It wraps the outputs of all required upstream
    stages for corresponding targets (e.g. all GVCFs from a GenotypeSample stage
    for a JointCalling stage, along with Hail Batch jobs).

    An object of this class is passed to the public `queue_jobs` method of a Stage,
    and can be used to query dependency files and jobs.
    """

    def __init__(self, stage: 'Stage'):
        self.stage = stage
        self._outputs_by_target_by_stage: dict[str, dict[str, StageOutput | None]] = {}

    def add_other_stage_output(self, output: StageOutput):
        """
        Add output from another stage run.
        """
        assert output.stage is not None, output
        if not output.target.active:
            return
        if not output.target.get_samples():
            return
        if not output.data and not output.jobs:
            return
        stage_name = output.stage.name
        target_id = output.target.target_id
        if stage_name not in self._outputs_by_target_by_stage:
            self._outputs_by_target_by_stage[stage_name] = dict()
        self._outputs_by_target_by_stage[stage_name][target_id] = output

    def _each(
        self,
        fun: Callable,
        stage: StageDecorator,
    ):
        if stage.__name__ not in [s.name for s in self.stage.required_stages]:
            raise PipelineError(
                f'{self.stage.name}: getting inputs from stage {stage.__name__}, '
                f'but {stage.__name__} is not listed in required_stages. '
                f'Consider adding it into the decorator: '
                f'@stage(required_stages=[{stage.__name__}])'
            )

        if stage.__name__ not in self._outputs_by_target_by_stage:
            raise PipelineError(
                f'No inputs from {stage.__name__} for {self.stage.name} found '
                'after skipping targets with missing inputs. '
                + (
                    'Check the logs if all samples were missing inputs from previous '
                    'stages, and consider changing `workflow/first_stage`'
                    if get_config()['workflow'].get('skip_samples_with_missing_input')
                    else ''
                )
            )

        return {
            trg: fun(result)
            for trg, result in self._outputs_by_target_by_stage.get(
                stage.__name__, {}
            ).items()
        }

    def as_path_by_target(
        self,
        stage: StageDecorator,
        id: str | None = None,
    ) -> dict[str, Path]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_path(id=id)), stage=stage)

    def as_resource_by_target(
        self,
        stage: StageDecorator,
        id: str | None = None,
    ) -> dict[str, hb.Resource]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_resource(id=id)), stage=stage)

    def as_dict_by_target(self, stage: StageDecorator) -> dict[str, dict[str, Path]]:
        """
        Get as a dict of files/resources for a specific stage, indexed by target
        """
        return self._each(fun=(lambda r: r.as_dict()), stage=stage)

    def as_resource_dict_by_target(
        self,
        stage: StageDecorator,
    ) -> dict[str, dict[str, hb.Resource]]:
        """
        Get a dict of resources for a specific stage, and indexed by target
        """
        return self._each(fun=(lambda r: r.as_resource_dict()), stage=stage)

    def as_path_dict_by_target(
        self,
        stage: StageDecorator,
    ) -> dict[str, dict[str, Path]]:
        """
        Get a dict of paths for a specific stage, and indexed by target
        """
        return self._each(fun=(lambda r: r.as_path_dict()), stage=stage)

    def _get(
        self,
        target: 'Target',
        stage: StageDecorator,
    ):
        if not self._outputs_by_target_by_stage.get(stage.__name__):
            raise StageInputNotFound(
                f'Not found output from stage {stage.__name__}, required for stage '
                f'{self.stage.name}'
            )
        if not self._outputs_by_target_by_stage[stage.__name__].get(target.target_id):
            raise StageInputNotFound(
                f'Not found output for {target} from stage {stage.__name__}, required '
                f'for stage {self.stage.name}'
            )
        return self._outputs_by_target_by_stage[stage.__name__][target.target_id]

    def as_path(
        self,
        target: 'Target',
        stage: StageDecorator,
        id: str | None = None,
    ) -> Path:
        """
        Represent as a path to a file, otherwise fail.
        `stage` can be callable, or a subclass of Stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_path(id)

    def as_resource(
        self,
        target: 'Target',
        stage: StageDecorator,
        id: str | None = None,
    ) -> hb.Resource:
        """
        Get Hail Batch Resource for a specific target and stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_resource(id)

    def as_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dict of files or Resources for a specific target and stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_dict()

    def as_path_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dict of files for a specific target and stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_path_dict()

    def as_resource_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, hb.Resource]:
        """
        Get a dict of  Resources for a specific target and stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_resource_dict()

    def get_jobs(self, target: 'Target') -> list[Job]:
        """
        Get list of jobs that the next stage would depend on.
        """
        all_jobs: list[Job] = []
        these_samples = target.get_sample_ids()
        for _, outputs_by_target in self._outputs_by_target_by_stage.items():
            for _, output in outputs_by_target.items():
                if output:
                    those_samples = output.target.get_sample_ids()
                    samples_intersect = set(these_samples) & set(those_samples)
                    if samples_intersect:
                        all_jobs.extend(output.jobs)
        return all_jobs


class Action(Enum):
    """
    Indicates what a stage should do with a specific target.
    """

    QUEUE = 1
    SKIP = 2
    REUSE = 3


class Stage(Generic[TargetT], ABC):
    """
    Abstract class for a pipeline stage. Parametrised by specific Target subclass,
    i.e. SampleStage(Stage[Sample]) should only be able to work on Sample(Target).
    """

    def __init__(
        self,
        name: str,
        batch: RegisteringBatch,
        cohort: Cohort,
        pipeline_tmp_prefix: Path,
        run_id: str,
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        analysis_type: str | None = None,
        skipped: bool = False,
        assume_outputs_exist: bool = False,
        forced: bool = False,
        status_reporter: StatusReporter | None = None,
    ):
        self._name = name
        self.b = batch
        self.cohort = cohort

        self.required_stages_classes: list[StageDecorator] = []
        if required_stages:
            if isinstance(required_stages, list):
                self.required_stages_classes.extend(required_stages)
            else:
                self.required_stages_classes.append(required_stages)

        self.tmp_prefix = pipeline_tmp_prefix / name
        self.run_id = run_id

        # Dependencies. Populated in pipeline.run(), after we know all stages.
        self.required_stages: list[Stage] = []

        self.status_reporter = status_reporter
        # If analysis type is defined, it will be used to update analysis status,
        # as well as find and reuse existing outputs from the status reporter
        self.analysis_type = analysis_type

        # Populated with the return value of `add_to_the_pipeline()`
        self.output_by_target: dict[str, StageOutput | None] = dict()

        self.skipped = skipped
        self.forced = forced
        self.assume_outputs_exist = assume_outputs_exist

    def __str__(self):
        res = f'{self._name}'
        if self.skipped:
            res += ' [skipped]'
        if self.forced:
            res += ' [forced]'
        if self.assume_outputs_exist:
            res += ' [assume_outputs_exist]'
        if self.required_stages:
            res += f' <- [{", ".join([s.name for s in self.required_stages])}]'
        return res

    @property
    def name(self):
        """
        Stage name (unique and descriptive stage)
        """
        return self._name

    @abstractmethod
    def queue_jobs(self, target: TargetT, inputs: StageInput) -> StageOutput | None:
        """
        Adds Hail Batch jobs that process `target`.
        Assumes that all the household work is done: checking missing inputs
        from required stages, checking for possible reuse of existing outputs.
        """

    @abstractmethod
    def expected_outputs(self, target: TargetT) -> ExpectedResultT:
        """
        Get path(s) to files that the stage is expected to generate for a `target`.
        Used within in `queue_jobs()` to pass paths to outputs to job commands,
        as well as by the pipeline to check if the stage's expected outputs already
        exist and can be reused.

        Can be a str, a Path object, or a dictionary of str/Path objects.
        """

    @abstractmethod
    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Queues jobs for each corresponding target, defined by Stage subclass.

        Returns a dictionary of `StageOutput` objects indexed by target unique_id.
        """

    def _make_inputs(self) -> StageInput:
        """
        Collects outputs from all dependencies and create input for this stage
        """
        inputs = StageInput(self)
        for prev_stage in self.required_stages:
            for _, stage_output in prev_stage.output_by_target.items():
                if stage_output:
                    inputs.add_other_stage_output(stage_output)
        return inputs

    def make_outputs(
        self,
        target: 'Target',
        data: StageOutputData | str | dict[str, str] | None = None,
        jobs: list[Job] | Job | None = None,
        meta: dict | None = None,
        reusable: bool = False,
        skipped: bool = False,
        error_msg: str | None = None,
    ) -> StageOutput:
        """
        Create StageOutput for this stage.
        """
        return StageOutput(
            target=target,
            data=data,
            jobs=jobs,
            meta=meta,
            reusable=reusable,
            skipped=skipped,
            error_msg=error_msg,
            stage=self,
        )

    def _queue_jobs_with_checks(
        self,
        target: TargetT,
        action: Action | None = None,
    ) -> StageOutput | None:
        """
        Checks what to do with target, and either queue jobs, or skip/reuse results.
        """
        if not action:
            action = self._get_action(target)

        inputs = self._make_inputs()
        expected_out = self.expected_outputs(target)

        if action == Action.QUEUE:
            outputs = self.queue_jobs(target, inputs)
        elif action == Action.REUSE:
            outputs = self.make_outputs(
                target=target,
                data=expected_out,
                reusable=True,
            )
        else:  # Action.SKIP
            outputs = None

        if not outputs:
            return None

        outputs.stage = self
        outputs.meta |= self.get_job_attrs(target)

        for j in outputs.jobs:
            if j:
                j.depends_on(*inputs.get_jobs(target))

        if outputs.error_msg:
            return outputs

        # Adding status reporter jobs
        if self.analysis_type and self.status_reporter:
            self.status_reporter.add_updaters_jobs(
                b=self.b,
                output=outputs.data,
                analysis_type=self.analysis_type,
                target=target,
                jobs=outputs.jobs,
                prev_jobs=inputs.get_jobs(target),
                meta=outputs.meta,
            )
        return outputs

    def new_reuse_job(self, target: Target) -> Job:
        """
        Add "reuse" job. Target doesn't have to be specific for a stage here,
        this using abstract class Target instead of generic parameter TargetT.
        """
        attrs = dict(stage=self.name, reuse=True)
        attrs |= target.get_job_attrs()
        return self.b.new_job(self.name, attrs)

    def _get_action(self, target: TargetT) -> Action:
        """
        Based on stage parameters and expected outputs existence, determines what
        to do with the target: queue, skip or reuse, etc..
        """
        if target.forced and not self.skipped:
            return Action.QUEUE

        if (
            d := get_config()['workflow'].get('skip_samples_stages')
        ) and self.name in d:
            skip_targets = d[self.name]
            if target.target_id in skip_targets:
                logger.info(f'{self.name}: requested to skip {target}')
                return Action.SKIP

        expected_out = self.expected_outputs(target)
        reusable, first_missing_path = self._is_reusable(expected_out)

        if self.skipped:
            if reusable and not first_missing_path:
                return Action.REUSE
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(
                    f'Skipping {target}: stage {self.name} is required, '
                    f'but is marked as skipped, and some expected outputs for the '
                    f'target do not exist: {first_missing_path}'
                )
                # `workflow/skip_samples_with_missing_input` means that we can ignore
                # samples/datasets that have missing results from skipped stages.
                # This is our case, so indicating that this sample/dataset should
                # be ignored:
                target.active = False
                return Action.SKIP
            if self.name in get_config()['workflow'].get(
                'allow_missing_outputs_for_stages', []
            ):
                return Action.REUSE
            else:
                raise ValueError(
                    f'Stage {self.name} is required, but is skipped, and '
                    f'the following expected outputs for target {target} do not exist: '
                    f'{first_missing_path}'
                )

        if reusable and not first_missing_path:
            if target.forced:
                logger.info(
                    f'{self.name}: can reuse, but forcing the target '
                    f'{target} to rerun this stage'
                )
                return Action.QUEUE
            elif self.forced:
                logger.info(
                    f'{self.name}: can reuse, but forcing the stage '
                    f'to rerun, target={target}'
                )
                return Action.QUEUE
            else:
                logger.info(f'{self.name}: reusing results for {target}')
                return Action.REUSE

        logger.info(f'{self.name}: running queue_jobs(target={target})')
        return Action.QUEUE

    def _is_reusable(self, expected_out: ExpectedResultT) -> tuple[bool, Path | None]:
        if self.assume_outputs_exist:
            return True, None

        if get_config()['workflow'].get('check_expected_outputs'):
            paths = []
            if isinstance(expected_out, Path):
                paths.append(expected_out)
            if isinstance(expected_out, dict):
                for _, v in expected_out.items():
                    if isinstance(v, Path):
                        paths.append(v)
            first_missing_path = next((p for p in paths if not exists(p)), None)
            if not paths:
                return False, None
            if first_missing_path:
                return False, first_missing_path
            return True, None
        else:
            if self.skipped:
                # Do not check the files' existence, trust they exist.
                # note that for skipped stages, we automatically assume outputs exist
                return True, None
            # Do not check the files' existence, assume they don't exist:
            return False, None

    def _queue_reuse_job(
        self, target: TargetT, found_paths: Path | dict[str, Path]
    ) -> StageOutput | None:
        """
        Queues a [reuse] Job
        """
        return self.make_outputs(
            target=target,
            data=found_paths,
            jobs=[self.b.new_job(f'{self.name} [reuse]', target.get_job_attrs())],
        )

    def get_job_attrs(self, target: TargetT | None = None) -> dict[str, str]:
        """
        Create Hail Batch Job attributes dictionary
        """
        job_attrs = dict(
            seq_type=get_config()['workflow']['sequencing_type'],
            stage=self.name,
        )
        if target:
            job_attrs |= target.get_job_attrs()
        return job_attrs


def stage(
    cls: Optional[Type['Stage']] = None,
    *,
    analysis_type: str | None = None,
    required_stages: list[StageDecorator] | StageDecorator | None = None,
    skipped: bool = False,
    assume_outputs_exist: bool = False,
    forced: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with optional arguments.
    The goal is to allow declaring pipeline stages without requiring to implement
    a constructor method. E.g.

    @stage(required_stages=[Align])
    class GenotypeSample(SampleStage):
        def expected_outputs(self, sample: Sample):
            ...
        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
            ...
    """

    def decorator_stage(_cls) -> StageDecorator:
        """Implements decorator."""

        @functools.wraps(_cls)
        def wrapper_stage(pipeline: 'Pipeline') -> Stage:
            """Decorator helper function."""
            return _cls(
                name=_cls.__name__,
                batch=pipeline.b,
                cohort=pipeline.cohort,
                pipeline_tmp_prefix=pipeline.tmp_prefix,
                run_id=pipeline.run_id,
                required_stages=required_stages,
                analysis_type=analysis_type,
                status_reporter=pipeline.status_reporter,
                skipped=skipped,
                assume_outputs_exist=assume_outputs_exist,
                forced=forced,
            )

        _ALL_DECLARED_STAGES.append(wrapper_stage)
        return wrapper_stage

    if cls is None:
        return decorator_stage
    else:
        return decorator_stage(cls)


# noinspection PyUnusedLocal
def skip(
    _fun: Optional[StageDecorator] = None,
    *,
    reason: str = None,
    assume_outputs_exist: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Decorator on top of `@stage` that sets the `self.skipped` field to True.
    By default, expected outputs of a skipped stage will be checked,
    unless `assume_outputs_exist` is True.

    @skip
    @stage
    class MyStage1(SampleStage):
        ...

    @skip
    @stage(assume_outputs_exist=True)
    class MyStage2(SampleStage):
        ...
    """

    def decorator_stage(fun) -> StageDecorator:
        """Implements decorator."""

        @functools.wraps(fun)
        def wrapper_stage(*args, **kwargs) -> Stage:
            """Decorator helper function."""
            s = fun(*args, **kwargs)
            s.skipped = True
            s.assume_outputs_exist = assume_outputs_exist
            return s

        return wrapper_stage

    if _fun is None:
        return decorator_stage
    else:
        return decorator_stage(_fun)


class Pipeline:
    """
    Represents a Pipeline, and encapsulates a Hail Batch object, stages,
    and a cohort of datasets of samples.
    """

    def __init__(
        self,
        name: str | None = None,
        description: str | None = None,
        stages: list[StageDecorator] | None = None,
    ):
        self._stages = stages
        self.run_id = get_config()['workflow'].get('run_id', timestamp())

        analysis_dataset = get_config()['workflow']['dataset']
        name = get_config()['workflow'].get('name') or name
        description = get_config()['workflow'].get('description') or description
        name = name or description or analysis_dataset
        self.name = slugify(name)
        description = description or name
        description += f': run_id={self.run_id}'

        self.cohort = get_cohort()
        self.tmp_prefix = self.cohort.analysis_dataset.tmp_prefix() / self.run_id
        description += f' [{get_config()["workflow"]["sequencing_type"]}]'
        if ds_set := set(d.name for d in self.cohort.get_datasets()):
            description += ' ' + ', '.join(sorted(ds_set))
        self.b: RegisteringBatch = setup_batch(description=description)

        self.status_reporter = None
        if get_config()['workflow'].get('status_reporter') == 'metamist':
            self.status_reporter = MetamistStatusReporter()

        # Will be populated by set_stages() in submit_batch()
        self._stages_dict: dict[str, Stage] = dict()

    def run(
        self,
        stages: list[StageDecorator] | None = None,
        wait: bool | None = False,
        force_all_implicit_stages: bool = False,
    ) -> Any | None:
        """
        Resolve stages, add and submit Hail Batch jobs.
        When `run_all_implicit_stages` is set, all required stages that were not defined
        explicitly would still be executed.
        """
        _stages_in_order = stages or self._stages or _ALL_DECLARED_STAGES
        if not _stages_in_order:
            raise PipelineError('No stages added')
        self.set_stages(_stages_in_order, force_all_implicit_stages)

        result = None
        if self.b:
            result = self.b.run(wait=wait)
        return result

    def _validate_first_last_stage(self) -> tuple[int | None, int | None]:
        """
        Validating the first_stage and the last_stage parameters.
        """
        first_stage = get_config()['workflow'].get('first_stage')
        last_stage = get_config()['workflow'].get('last_stage')

        stage_names = list(self._stages_dict.keys())
        lower_stage_names = [s.lower() for s in stage_names]
        first_stage_num = None
        if first_stage:
            if first_stage.lower() not in lower_stage_names:
                logger.critical(
                    f'Value for --first-stage {first_stage} '
                    f'not found in available stages: {", ".join(stage_names)}'
                )
            first_stage_num = lower_stage_names.index(first_stage.lower())
        last_stage_num = None
        if last_stage:
            if last_stage.lower() not in lower_stage_names:
                logger.critical(
                    f'Value for --last-stage {last_stage} '
                    f'not found in available stages: {", ".join(stage_names)}'
                )
            last_stage_num = lower_stage_names.index(last_stage.lower())
        return first_stage_num, last_stage_num

    def set_stages(
        self,
        stages_classes: list[StageDecorator],
        force_all_implicit_stages: bool = False,
    ):
        """
        Iterate over stages and call queue_for_cohort(cohort) on each.
        Effectively creates all Hail Batch jobs through Stage.queue_jobs().

        When `run_all_implicit_stages` is set, all required stages that were not set
        explicitly would still be run.
        """
        if not self.cohort:
            raise PipelineError('Cohort must be populated before adding stages')

        # First round: initialising stage objects
        stages = [cls(self) for cls in stages_classes]
        for stage_ in stages:
            if stage_.name in self._stages_dict:
                raise ValueError(
                    f'Stage {stage_.name} is already defined. Check this '
                    f'list for duplicates: {", ".join(s.name for s in stages)}'
                )
            self._stages_dict[stage_.name] = stage_

        # Second round: checking which stages are required, even implicitly.
        # implicit_stages_d: dict[str, Stage] = dict()  # If there are required
        # dependency stages that are not requested explicitly, we are putting them
        # into this dict as `skipped`.
        depth = 0
        while True:  # Might need several rounds to resolve dependencies recursively.
            depth += 1
            newly_implicitly_added_d = dict()
            for stage_ in self._stages_dict.values():
                for reqcls in stage_.required_stages_classes:  # check dependencies
                    if reqcls.__name__ in self._stages_dict:  # already added
                        continue
                    # Initialising and adding as explicit.
                    reqstage = reqcls(self)
                    newly_implicitly_added_d[reqstage.name] = reqstage
                    if reqcls.__name__ in get_config()['workflow'].get(
                        'assume_outputs_exist_for_stages', []
                    ):
                        reqstage.assume_outputs_exist = True
                    if reqcls.__name__ in get_config()['workflow'].get(
                        'skip_stages', []
                    ):
                        reqstage.skipped = True
                        continue
                    if not force_all_implicit_stages:
                        # Stage is not declared or requested implicitly, so setting
                        # it as skipped:
                        reqstage.skipped = True
                        # Only checking outputs of immediately required stages
                        if depth > 1:
                            reqstage.assume_outputs_exist = True
                            logger.info(f'Stage {reqstage.name} is skipped')
                        else:
                            logger.info(
                                f'Stage {reqstage.name} is skipped, but the output '
                                f'will be required for the stage {stage_.name}'
                            )

            if newly_implicitly_added_d:
                logger.info(
                    f'Additional implicit stages: '
                    f'{list(newly_implicitly_added_d.keys())}'
                )
                # Adding new stages back into the ordered dict, so they are
                # executed first.
                self._stages_dict = newly_implicitly_added_d | self._stages_dict
            else:
                # No new implicit stages added, can stop here.
                break

        for stage_ in self._stages_dict.values():
            stage_.required_stages = [
                self._stages_dict[cls.__name__]
                for cls in stage_.required_stages_classes
            ]

        # Second round - applying first and last stage options.
        first_stage_num, last_stage_num = self._validate_first_last_stage()
        for i, (stage_name, stage_) in enumerate(self._stages_dict.items()):
            if first_stage_num is not None and i < first_stage_num:
                stage_.skipped = True
                if i < first_stage_num - 1:
                    # Not checking expected outputs of stages before that
                    stage_.assume_outputs_exist = True
                logger.info(f'Skipping stage {stage_name}')
                continue
            if last_stage_num is not None and i > last_stage_num:
                stage_.skipped = True
                stage_.assume_outputs_exist = True
                continue

        if not (
            final_set_of_stages := [
                s.name for s in self._stages_dict.values() if not s.skipped
            ]
        ):
            raise PipelineError('No stages to run')
        logger.info(f'Setting stages: {final_set_of_stages}')
        required_skipped_stages = [s for s in self._stages_dict.values() if s.skipped]
        if required_skipped_stages:
            logger.info(
                f'Skipped stages: ' f'{[s.name for s in required_skipped_stages]}'
            )

        # Second round - actually adding jobs from the stages.
        for i, (_, stage_) in enumerate(self._stages_dict.items()):
            logger.info(f'*' * 60)
            logger.info(f'Stage {stage_}')
            stage_.output_by_target = stage_.queue_for_cohort(self.cohort)
            if errors := self._process_stage_errors(stage_.output_by_target):
                raise PipelineError(
                    f'Stage {stage_} failed to queue jobs with errors: '
                    + '\n'.join(errors)
                )

            logger.info(f'')
            if last_stage_num is not None and i >= last_stage_num:
                logger.info(f'Last stage was {stage_.name}, stopping here')
                break

    @staticmethod
    def _process_stage_errors(
        output_by_target: dict[str, StageOutput | None]
    ) -> list[str]:
        targets_by_error = defaultdict(list)
        for target, output in output_by_target.items():
            if output and output.error_msg:
                targets_by_error[output.error_msg].append(target)
        return [
            f'{error}: {", ".join(target_ids)}'
            for error, target_ids in targets_by_error.items()
        ]

    def get_datasets(self, only_active: bool = True) -> list[Dataset]:
        """
        Thin wrapper.
        """
        return self.cohort.get_datasets(only_active=only_active)

    def create_dataset(self, name: str) -> Dataset:
        """
        Thin wrapper.
        """
        return self.cohort.create_dataset(name)

    def get_all_samples(self, only_active: bool = True) -> list[Sample]:
        """
        Thin wrapper.
        """
        return self.cohort.get_samples(only_active=only_active)

    def get_all_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Thin wrapper.
        """
        return self.cohort.get_sample_ids(only_active=only_active)
