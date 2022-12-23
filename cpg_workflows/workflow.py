"""
Provides a `Workflow` class and a `@stage` decorator that allow to define workflows
in a declarative fashion.

A `Stage` object is responsible for creating Hail Batch jobs and declaring outputs
(files or metamist analysis objects) that are expected to be produced. Each stage
acts on a `Target`, which can be of the following a `Sample`, a `Dataset` (a container
for samples), or a `Cohort` (all input datasets together). A `Workflow` object plugs
stages together by resolving dependencies between different levels accordingly.

Examples of workflows can be found in the `production-workflows` repository.
"""

import functools
import networkx as nx
import logging
import pathlib
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from typing import cast, Callable, Union, TypeVar, Generic, Optional, Type, Sequence

from cloudpathlib import CloudPath
from hailtop.batch.job import Job
from cpg_utils.config import get_config
from cpg_utils import Path

from .batch import get_batch
from .status import MetamistStatusReporter
from .targets import Target, Dataset, Sample, Cohort
from .utils import exists, timestamp, slugify
from .inputs import get_cohort


StageDecorator = Callable[..., 'Stage']

# Type variable to use with Generic to make sure a Stage subclass always matches the
# corresponding Target subclass. We can't just use the Target superclass because
# it would violate the Liskov substitution principle (i.e. any Stage subclass would
# have to be able to work on any Target subclass).
TargetT = TypeVar('TargetT', bound=Target)

ExpectedResultT = Union[Path, dict[str, Path], dict[str, str], str, None]

StageOutputData = Union[Path, dict[str, Path]]


class WorkflowError(Exception):
    """
    Error raised by workflow and stage implementation.
    """


class StageInputNotFoundError(Exception):
    """
    Thrown when a stage requests input from another stage
    that doesn't exist.
    """


# noinspection PyShadowingNames
class StageOutput:
    """
    Represents a result of a specific stage, which was run on a specific target.
    Can be a file path, or a Hail Batch Resource. Optionally wrapped in a dict.
    """

    def __init__(
        self,
        target: 'Target',
        data: StageOutputData | None = None,
        jobs: Sequence[Job | None] | Job | None = None,
        meta: dict | None = None,
        reusable: bool = False,
        skipped: bool = False,
        error_msg: str | None = None,
        stage: Optional['Stage'] = None,
    ):
        # Converting str into Path objects.
        self.data = data
        self.stage = stage
        self.target = target
        _jobs = [jobs] if isinstance(jobs, Job) else (jobs or [])
        self.jobs: list[Job] = [j for j in _jobs if j is not None]
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

    def _get(self, key=None) -> str | Path:
        if self.data is None:
            raise ValueError(f'{self.stage}: output data is not available')

        if key is not None:
            if not isinstance(self.data, dict):
                raise ValueError(
                    f'{self.stage}: {self.data} is not a dictionary, can\'t get "{key}"'
                )
            res = cast(dict, self.data)[key]
        else:
            res = self.data
        return res

    def as_str(self, key=None) -> str:
        """
        Cast the result to a simple string. Throw an exception when can't cast.
        `key` is used to extract the value when the result is a dictionary.
        """
        res = self._get(key)
        if not isinstance(res, str):
            raise ValueError(f'{res} is not a str.')
        return cast(str, res)

    def as_path(self, key=None) -> Path:
        """
        Cast the result to a path object. Throw an exception when can't cast.
        `key` is used to extract the value when the result is a dictionary.
        """
        res = self._get(key)
        if not isinstance(res, CloudPath | pathlib.Path):
            raise ValueError(f'{res} is not a path object.')

        return cast(Path, res)

    def as_dict(self) -> dict[str, Path]:
        """
        Cast the result to a dictionary, or throw an error if the cast failed.
        """
        if not isinstance(self.data, dict):
            raise ValueError(f'{self.data} is not a dictionary.')
        return self.data


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
            raise WorkflowError(
                f'{self.stage.name}: getting inputs from stage {stage.__name__}, '
                f'but {stage.__name__} is not listed in required_stages. '
                f'Consider adding it into the decorator: '
                f'@stage(required_stages=[{stage.__name__}])'
            )

        if stage.__name__ not in self._outputs_by_target_by_stage:
            raise WorkflowError(
                f'No inputs from {stage.__name__} for {self.stage.name} found '
                + 'after skipping targets with missing inputs. '
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
        key: str | None = None,
    ) -> dict[str, Path]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_path(key=key)), stage=stage)

    def as_dict_by_target(self, stage: StageDecorator) -> dict[str, dict[str, Path]]:
        """
        Get as a dict of files/resources for a specific stage, indexed by target
        """
        return self._each(fun=(lambda r: r.as_dict()), stage=stage)

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
            raise StageInputNotFoundError(
                f'Not found output from stage {stage.__name__}, required for stage '
                f'{self.stage.name}. Available: {self._outputs_by_target_by_stage}'
            )
        if not self._outputs_by_target_by_stage[stage.__name__].get(target.target_id):
            raise StageInputNotFoundError(
                f'Not found output for {target} from stage {stage.__name__}, required '
                f'for stage {self.stage.name}'
            )
        return self._outputs_by_target_by_stage[stage.__name__][target.target_id]

    def as_path(
        self,
        target: 'Target',
        stage: StageDecorator,
        key: str | None = None,
    ) -> Path:
        """
        Represent as a path to a file, otherwise fail.
        `stage` can be callable, or a subclass of Stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_path(key)

    def as_str(
        self,
        target: 'Target',
        stage: StageDecorator,
        key: str | None = None,
    ) -> str:
        """
        Represent as a simple string, otherwise fail.
        `stage` can be callable, or a subclass of Stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_str(key)

    def as_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dict of paths for a specific target and stage
        """
        res = self._get(target=target, stage=stage)
        return res.as_dict()

    def get_jobs(self, target: 'Target') -> list[Job]:
        """
        Get list of jobs that the next stage would depend on.
        """
        all_jobs: list[Job] = []
        these_samples = target.get_sample_ids()
        for stage_, outputs_by_target in self._outputs_by_target_by_stage.items():
            for target_, output in outputs_by_target.items():
                if output:
                    those_samples = output.target.get_sample_ids()
                    samples_intersect = set(these_samples) & set(those_samples)
                    if samples_intersect:
                        for j in output.jobs:
                            assert (
                                j
                            ), f'Stage: {stage_}, target: {target_}, output: {output}'
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
    Abstract class for a workflow stage. Parametrised by specific Target subclass,
    i.e. SampleStage(Stage[Sample]) should only be able to work on Sample(Target).
    """

    def __init__(
        self,
        name: str,
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        analysis_type: str | None = None,
        analysis_key: str | None = None,
        skipped: bool = False,
        assume_outputs_exist: bool = False,
        forced: bool = False,
    ):
        self._name = name
        self.required_stages_classes: list[StageDecorator] = []
        if required_stages:
            if isinstance(required_stages, list):
                self.required_stages_classes.extend(required_stages)
            else:
                self.required_stages_classes.append(required_stages)

        # Dependencies. Populated in workflow.run(), after we know all stages.
        self.required_stages: list[Stage] = []

        self.status_reporter = get_workflow().status_reporter
        # If `analysis_type` is defined, it will be used to create/update Analysis
        # entries in Metamist.
        self.analysis_type = analysis_type
        # If `analysis_key` is defined, it will be used to extract the value for
        # `Analysis.output` if the Stage.expected_outputs() returns a dict.
        self.analysis_key = analysis_key

        # Populated with the return value of `add_to_the_workflow()`
        self.output_by_target: dict[str, StageOutput | None] = dict()

        self.skipped = skipped
        self.forced = forced or self.name in get_config()['workflow'].get(
            'force_stages', []
        )
        self.assume_outputs_exist = assume_outputs_exist

    @property
    def tmp_prefix(self):
        return get_workflow().tmp_prefix / self.name

    @property
    def web_prefix(self) -> Path:
        return get_workflow().web_prefix / self.name

    @property
    def prefix(self) -> Path:
        return get_workflow().prefix / self.name

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
    def name(self) -> str:
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
        as well as by the workflow to check if the stage's expected outputs already
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
        jobs: Sequence[Job | None] | Job | None = None,
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

        for output_job in outputs.jobs:
            if output_job:
                for input_job in inputs.get_jobs(target):
                    assert (
                        input_job
                    ), f'Input dependency job for stage: {self}, target: {target}'
                    output_job.depends_on(input_job)

        if outputs.error_msg:
            return outputs

        # Adding status reporter jobs
        if (
            self.analysis_type
            and self.status_reporter
            and action == Action.QUEUE
            and outputs.data
        ):
            output: str | Path
            if isinstance(outputs.data, dict):
                if not self.analysis_key:
                    raise WorkflowError(
                        f'Cannot create Analysis for stage {self.name}: `analysis_key` '
                        f'must be set with the @stage decorator to select value from '
                        f'the expected_outputs dict: {outputs.data}'
                    )
                if self.analysis_key not in outputs.data:
                    raise WorkflowError(
                        f'Cannot create Analysis for stage {self.name}: `analysis_key` '
                        f'"{self.analysis_key}" is not found in the expected_outputs '
                        f'dict {outputs.data}'
                    )
                output = outputs.data[self.analysis_key]
            else:
                output = outputs.data

            assert isinstance(output, str) or isinstance(output, Path), output

            self.status_reporter.add_updaters_jobs(
                b=get_batch(),
                output=str(output),
                analysis_type=self.analysis_type,
                target=target,
                jobs=outputs.jobs,
                prev_jobs=inputs.get_jobs(target),
                meta=outputs.meta,
                job_attrs=self.get_job_attrs(target),
            )
        return outputs

    def _get_action(self, target: TargetT) -> Action:
        """
        Based on stage parameters and expected outputs existence, determines what
        to do with the target: queue, skip or reuse, etc..
        """
        if target.forced and not self.skipped:
            logging.info(f'{self.name}: {target} [QUEUE] (target is forced)')
            return Action.QUEUE

        if (
            d := get_config()['workflow'].get('skip_samples_stages')
        ) and self.name in d:
            skip_targets = d[self.name]
            if target.target_id in skip_targets:
                logging.info(
                    f'{self.name}: {target} [SKIP] (is in workflow/skip_samples_stages)'
                )
                return Action.SKIP

        expected_out = self.expected_outputs(target)
        reusable, first_missing_path = self._is_reusable(expected_out)

        if self.skipped:
            if reusable and not first_missing_path:
                logging.info(
                    f'{self.name}: {target} [REUSE] (stage skipped, and outputs exist)'
                )
                return Action.REUSE
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logging.warning(
                    f'{self.name}: {target} [SKIP] (stage is required, '
                    f'but is marked as "skipped", '
                    f'workflow/skip_samples_with_missing_input=true '
                    f'and some expected outputs for the target do not exist: '
                    f'{first_missing_path}'
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
                logging.info(
                    f'{self.name}: {target} [REUSE] (stage is skipped, some outputs are'
                    f'missing, but stage is listed in '
                    f'workflow/allow_missing_outputs_for_stages)'
                )
                return Action.REUSE
            else:
                raise WorkflowError(
                    f'{self.name}: stage is required, but is skipped, and '
                    f'the following expected outputs for target {target} do not exist: '
                    f'{first_missing_path}'
                )

        if reusable and not first_missing_path:
            if target.forced:
                logging.info(
                    f'{self.name}: {target} [QUEUE] (can reuse, but forcing the target '
                    f'to rerun this stage)'
                )
                return Action.QUEUE
            elif self.forced:
                logging.info(
                    f'{self.name}: {target} [QUEUE] (can reuse, but forcing the stage '
                    f'to rerun)'
                )
                return Action.QUEUE
            else:
                logging.info(
                    f'{self.name}: {target} [REUSE] (expected outputs exist: '
                    f'{expected_out})'
                )
                return Action.REUSE

        logging.info(f'{self.name}: {target} [QUEUE]')
        return Action.QUEUE

    def _is_reusable(self, expected_out: ExpectedResultT) -> tuple[bool, Path | None]:
        if self.assume_outputs_exist:
            return True, None

        if not expected_out:
            # Marking is reusable. If the stage does not naturally produce any outputs,
            # it would still need to create some flag file.
            return True, None

        if get_config()['workflow'].get('check_expected_outputs'):
            paths = []
            if isinstance(expected_out, dict):
                for _, v in expected_out.items():
                    if not isinstance(v, str):
                        paths.append(v)
            elif isinstance(expected_out, str):
                pass
            else:
                paths.append(expected_out)
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

    def get_job_attrs(self, target: TargetT | None = None) -> dict[str, str]:
        """
        Create Hail Batch Job attributes dictionary
        """
        job_attrs = dict(stage=self.name)
        if sequencing_type := get_config()['workflow'].get('sequencing_type'):
            job_attrs['sequencing_type'] = sequencing_type
        if target:
            job_attrs |= target.get_job_attrs()
        return job_attrs


def stage(
    cls: Optional[Type['Stage']] = None,
    *,
    analysis_type: str | None = None,
    analysis_key: str | None = None,
    required_stages: list[StageDecorator] | StageDecorator | None = None,
    skipped: bool = False,
    assume_outputs_exist: bool = False,
    forced: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with optional arguments.
    The goal is to allow declaring workflow stages without requiring to implement
    a constructor method. E.g.

    @stage(required_stages=[Align])
    class GenotypeSample(SampleStage):
        def expected_outputs(self, sample: Sample):
            ...
        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
            ...

        self.status_reporter = get_workflow().status_reporter
        # If `analysis_type` is defined, it will be used to create/update Analysis
        # entries in Metamist.
        self.analysis_type = analysis_type
        # If `analysis_key` is defined, it will be used to extract the value for
        # `Analysis.output` if the Stage.expected_outputs() returns a dict.

    @analysis_type: if defined, will be used to create/update `Analysis` entries
        using the status reporter.
    @analysis_key: is defined, will be used to extract the value for `Analysis.output`
        if the Stage.expected_outputs() returns a dict.
    @required_stages: list of other stage classes that are required prerequisites
        for this stage. Outputs of those stages will be passed to
        `Stage.queue_jobs(... , inputs)` as `inputs`, and all required
        dependencies between Hail Batch jobs will be set automatically as well.
    @skipped: always skip this stage.
    @assume_outputs_exist: assume expected outputs of this stage always exist.
    @forced: always force run this stage, regardless of the outputs' existence.
    """

    def decorator_stage(_cls) -> StageDecorator:
        """Implements decorator."""

        @functools.wraps(_cls)
        def wrapper_stage() -> Stage:
            """Decorator helper function."""
            return _cls(
                name=_cls.__name__,
                required_stages=required_stages,
                analysis_type=analysis_type,
                analysis_key=analysis_key,
                skipped=skipped,
                assume_outputs_exist=assume_outputs_exist,
                forced=forced,
            )

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


_workflow: Optional['Workflow'] = None


def get_workflow(dry_run: bool = False) -> 'Workflow':
    global _workflow
    if _workflow is None:
        _workflow = Workflow(dry_run=dry_run)
    return _workflow


def run_workflow(
    stages: list[StageDecorator] | None = None,
    wait: bool | None = False,
    dry_run: bool = False,
) -> 'Workflow':
    wfl = get_workflow(dry_run=dry_run)
    wfl = wfl.run(stages=stages, wait=wait)
    return wfl


class Workflow:
    """
    Encapsulates a Hail Batch object, stages, and a cohort of datasets of samples.
    Responsible for orchestrating stages.
    """

    def __init__(
        self,
        stages: list[StageDecorator] | None = None,
        dry_run: bool | None = None,
    ):
        if _workflow is not None:
            raise ValueError(
                'Workflow already initialised. Use get_workflow() to get the instance'
            )

        self.dry_run = dry_run or get_config()['workflow'].get('dry_run')

        analysis_dataset = get_config()['workflow']['dataset']
        name = get_config()['workflow'].get('name')
        description = get_config()['workflow'].get('description')
        name = name or description or analysis_dataset
        self.name = slugify(name)

        self._output_version: str | None = None
        if output_version := get_config()['workflow'].get('output_version'):
            self._output_version = slugify(output_version)

        self.run_timestamp: str = (
            get_config()['workflow'].get('run_timestamp') or timestamp()
        )

        # Description
        description = description or name
        if self._output_version:
            description += f': output_version={self._output_version}'
        description += f': run_timestamp={self.run_timestamp}'
        if sequencing_type := get_config()['workflow'].get('sequencing_type'):
            description += f' [{sequencing_type}]'
        if not self.dry_run:
            if ds_set := set(d.name for d in get_cohort().get_datasets()):
                description += ' ' + ', '.join(sorted(ds_set))
            get_batch().name = description

        self.status_reporter = None
        if get_config()['workflow'].get('status_reporter') == 'metamist':
            self.status_reporter = MetamistStatusReporter()
        self._stages: list[StageDecorator] | None = stages

    @property
    def output_version(self) -> str:
        return self._output_version or get_cohort().alignment_inputs_hash()

    @property
    def tmp_prefix(self) -> Path:
        return self._prefix(category='tmp')

    @property
    def web_prefix(self) -> Path:
        return self._prefix(category='web')

    @property
    def prefix(self) -> Path:
        return self._prefix()

    def _prefix(self, category: str | None = None) -> Path:
        """
        Prepare a unique path for the workflow with this name and this input data.
        """
        prefix = get_cohort().analysis_dataset.prefix(category=category) / self.name
        prefix /= self.output_version
        return prefix

    def run(
        self,
        stages: list[StageDecorator] | None = None,
        wait: bool | None = False,
    ):
        """
        Resolve stages, add and submit Hail Batch jobs.
        When `run_all_implicit_stages` is set, all required stages that were not defined
        explicitly would still be executed.
        """
        _stages = stages or self._stages
        if not _stages:
            raise WorkflowError('No stages added')
        self.set_stages(_stages)

        if not self.dry_run:
            get_batch().run(wait=wait)

    @staticmethod
    def _process_first_last_stages(
        stages: list[Stage],
        graph: nx.DiGraph,
        first_stages: list[str],
        last_stages: list[str],
    ):
        """
        Applying first_stages and last_stages config options. Would skip all stages
        before first_stages, and all stages after last_stages (i.e. descendants and
        ancestors on the stages DAG.)
        """
        stages_d = {s.name: s for s in stages}
        stage_names = list(stg.name for stg in stages)
        lower_names = {s.lower() for s in stage_names}

        for param, _stage_list in [
            ('first_stages', first_stages),
            ('last_stages', last_stages),
        ]:
            for _s_name in _stage_list:
                if _s_name.lower() not in lower_names:
                    raise WorkflowError(
                        f'Value in workflow/{param} "{_s_name}" must be a stage name '
                        f'or a subset of stages from the available list: '
                        f'{", ".join(stage_names)}'
                    )

        if not (last_stages or first_stages):
            return

        # E.g. if our last_stages is CramQc, MtToEs would still run because it's in
        # a different branch. So we want to collect all stages after first_stages
        # and before last_stages in their respective branches, and mark as skipped
        # everything in other branches.
        first_stages_keeps: list[str] = first_stages[:]
        last_stages_keeps: list[str] = last_stages[:]

        for fs in first_stages:
            for descendant in nx.descendants(graph, fs):
                if not stages_d[descendant].skipped:
                    logging.info(
                        f'Skipping stage {descendant} (precedes {fs} listed in '
                        f'first_stages)'
                    )
                    stages_d[descendant].skipped = True
                for grand_descendant in nx.descendants(graph, descendant):
                    if not stages_d[grand_descendant].assume_outputs_exist:
                        logging.info(
                            f'Not checking expected outputs of not immediately '
                            f'required stage {grand_descendant} (< {descendant} < {fs})'
                        )
                        stages_d[grand_descendant].assume_outputs_exist = True

            for ancestor in nx.ancestors(graph, fs):
                first_stages_keeps.append(ancestor)

        for ls in last_stages:
            for ancestor in nx.ancestors(graph, ls):
                if not stages_d[ancestor].skipped:
                    logging.info(f'Skipping stage {ancestor} (after last {ls})')
                    stages_d[ancestor].skipped = True
                    stages_d[ancestor].assume_outputs_exist = True

            for ancestor in nx.descendants(graph, ls):
                last_stages_keeps.append(ancestor)

        for _stage in stages:
            if _stage.name not in last_stages_keeps + first_stages_keeps:
                _stage.skipped = True
                _stage.assume_outputs_exist = True

    def set_stages(
        self,
        requested_stages: list[StageDecorator],
    ):
        """
        Iterate over stages and call their queue_for_cohort(cohort) methods;
        through that, creates all Hail Batch jobs through Stage.queue_jobs().
        """
        # TOML options to configure stages:
        skip_stages = get_config()['workflow'].get('skip_stages', [])
        only_stages = get_config()['workflow'].get('only_stages', [])
        first_stages = get_config()['workflow'].get('first_stages', [])
        last_stages = get_config()['workflow'].get('last_stages', [])

        logging.info(
            f'End stages for the workflow "{self.name}": '
            f'{[cls.__name__ for cls in requested_stages]}'
        )
        logging.info('Stages additional configuration:')
        logging.info(f'  workflow/skip_stages: {skip_stages}')
        logging.info(f'  workflow/only_stages: {only_stages}')
        logging.info(f'  workflow/first_stages: {first_stages}')
        logging.info(f'  workflow/last_stages: {last_stages}')

        # Round 1: initialising stage objects.
        _stages_d: dict[str, Stage] = {}
        for cls in requested_stages:
            if cls.__name__ in _stages_d:
                continue
            _stages_d[cls.__name__] = cls()

        # Round 2: depth search to find implicit stages.
        depth = 0
        while True:  # might require few iterations to resolve dependencies recursively
            depth += 1
            newly_implicitly_added_d = dict()
            for stg in _stages_d.values():
                if stg.name in skip_stages:
                    stg.skipped = True
                    continue  # not searching deeper
                if only_stages and stg.name not in only_stages:
                    stg.skipped = True

                # Iterate dependencies:
                for reqcls in stg.required_stages_classes:
                    if reqcls.__name__ in _stages_d:  # already added
                        continue
                    # Initialising and adding as explicit.
                    reqstg = reqcls()
                    newly_implicitly_added_d[reqstg.name] = reqstg

            if newly_implicitly_added_d:
                logging.info(
                    f'Additional implicit stages: '
                    f'{list(newly_implicitly_added_d.keys())}'
                )
                _stages_d |= newly_implicitly_added_d
            else:
                # No new implicit stages added, so can stop the depth-search here
                break

        # Round 3: set "stage.required_stages" fields to each stage.
        for stg in _stages_d.values():
            stg.required_stages = [
                _stages_d[cls.__name__]
                for cls in stg.required_stages_classes
                if cls.__name__ in _stages_d
            ]

        # Round 4: determining order of execution.
        dag_node2nodes = dict()  # building a DAG
        for stg in _stages_d.values():
            dag_node2nodes[stg.name] = set(dep.name for dep in stg.required_stages)
        dag = nx.DiGraph(dag_node2nodes)
        try:
            stage_names = list(reversed(list(nx.topological_sort(dag))))
        except nx.NetworkXUnfeasible:
            logging.error('Circular dependencies found between stages')
            raise
        logging.info(f'Stages in order of execution:\n{stage_names}')
        stages = [_stages_d[name] for name in stage_names]

        # Round 5: applying workflow options first_stages and last_stages.
        self._process_first_last_stages(stages, dag, first_stages, last_stages)

        if not (final_set_of_stages := [s.name for s in stages if not s.skipped]):
            raise WorkflowError('No stages to run')
        logging.info(
            f'Final list of stages after applying workflow/first_stages and '
            f'workflow/last_stages stages:\n{final_set_of_stages}'
        )
        required_skipped_stages = [s for s in stages if s.skipped]
        if required_skipped_stages:
            logging.info(
                f'Skipped stages: {", ".join(s.name for s in required_skipped_stages)}'
            )

        # Round 6: actually adding jobs from the stages.
        if not self.dry_run:
            cohort = get_cohort()  # Would communicate with metamist.
            for i, stg in enumerate(stages):
                logging.info(f'*' * 60)
                logging.info(f'Stage #{i + 1}: {stg}')
                stg.output_by_target = stg.queue_for_cohort(cohort)
                if errors := self._process_stage_errors(stg.output_by_target):
                    raise WorkflowError(
                        f'Stage {stg} failed to queue jobs with errors: '
                        + '\n'.join(errors)
                    )

                logging.info(f'')

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


class SampleStage(Stage[Sample], ABC):
    """
    Sample-level stage.
    """

    @abstractmethod
    def expected_outputs(self, sample: Sample) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()

        if not (datasets := cohort.get_datasets()):
            logging.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.input_datasets` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`'
            )
            return output_by_target
        if not cohort.get_samples():
            logging.warning(
                f'{len(cohort.get_samples())}/'
                f'{len(cohort.get_samples(only_active=False))} '
                f'usable (active=True) samples found. Check logs above for '
                f'possible reasons samples were skipped (e.g. all samples ignored '
                f'via `workflow.skip_samples` in config, or they all missing stage '
                f'inputs and `workflow.skip_samples_with_missing_input=true` is set)'
            )
            return output_by_target

        for dataset in datasets:
            if not dataset.get_samples():
                logging.warning(
                    f'{dataset}: '
                    f'{len(dataset.get_samples())}/'
                    f'{len(dataset.get_samples(only_active=False))} '
                    f'usable (active=True) samples found. Check logs above for '
                    f'possible reasons samples were skipped (e.g. all samples ignored '
                    f'via `workflow.skip_samples` in config, or they all missing stage '
                    f'inputs and `workflow.skip_samples_with_missing_input=true` is set)'
                )
                continue

            logging.info(f'Dataset {dataset}:')
            for sample in dataset.get_samples():
                action = self._get_action(sample)
                output_by_target[sample.target_id] = self._queue_jobs_with_checks(
                    sample, action
                )

        return output_by_target


class DatasetStage(Stage, ABC):
    """
    Dataset-level stage
    """

    @abstractmethod
    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        if not (datasets := cohort.get_datasets()):
            logging.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.input_datasets` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`'
            )
            return output_by_target
        for dataset_i, dataset in enumerate(datasets):
            action = self._get_action(dataset)
            logging.info(f'{self.name}: #{dataset_i + 1}/{dataset} [{action.name}]')
            output_by_target[dataset.target_id] = self._queue_jobs_with_checks(
                dataset, action
            )
        return output_by_target


class CohortStage(Stage, ABC):
    """
    Cohort-level stage (all datasets of a workflow run).
    """

    @abstractmethod
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        return {cohort.target_id: self._queue_jobs_with_checks(cohort)}
