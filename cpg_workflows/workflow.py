"""
Provides a `Workflow` class and a `@stage` decorator that allow to define workflows
in a declarative fashion.

A `Stage` object is responsible for creating Hail Batch jobs and declaring outputs
(files or metamist analysis objects) that are expected to be produced. Each stage
acts on a `Target`, which can be of the following a `SequencingGroup`, a `Dataset` (a container
for sequencing groups), or a `Cohort` (all input datasets together). A `Workflow` object plugs
stages together by resolving dependencies between different levels accordingly.

Examples of workflows can be found in the `production-workflows` repository.
"""

import functools
import logging
import pathlib
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from typing import Callable, Generic, Optional, Sequence, Type, TypeVar, Union, cast

import networkx as nx
from cloudpathlib import CloudPath

from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, reset_batch

from .inputs import get_multicohort
from .status import MetamistStatusReporter
from .targets import Cohort, Dataset, MultiCohort, SequencingGroup, Target
from .utils import (
    ExpectedResultT,
    exists,
    slugify,
    timestamp,
)

StageDecorator = Callable[..., 'Stage']

# Type variable to use with Generic to make sure a Stage subclass always matches the
# corresponding Target subclass. We can't just use the Target superclass because
# it would violate the Liskov substitution principle (i.e. any Stage subclass would
# have to be able to work on any Target subclass).
TargetT = TypeVar('TargetT', bound=Target)


def path_walk(expected, collected: set | None = None) -> set[Path]:
    """
    recursive walk of expected_out
    if the object is iterable, walk it
    this gets around the issue with nested lists and dicts
    mainly around the use of Array outputs from Cromwell

    Args:
        expected (Any): any type of object containing Paths
        collected (set): all collected paths so far

    Returns:
        a set of all collected Path nodes

    Examples:

    >>> path_walk({'a': {'b': {'c': Path('d')}}})
    {Path('d')}
    >>> path_walk({'a': {'b': {'c': [Path('d'), Path('e')]}}})
    {Path('d'), Path('e')}
    >>> path_walk({'a': Path('b'),'c': {'d': 'e'}, {'f': Path('g')}})
    {Path('b'), Path('g')}
    """
    if collected is None:
        collected = set()

    if expected is None:
        return collected
    if isinstance(expected, dict):
        for value in expected.values():
            collected.update(path_walk(value, collected))
    if isinstance(expected, (list, set)):
        for value in expected:
            collected.update(path_walk(value, collected))
    if isinstance(expected, str):
        return collected
    if isinstance(expected, Path):
        if expected in collected:
            raise ValueError(f'Duplicate path {expected} in expected_out')
        collected.add(expected)
    return collected


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
        data: ExpectedResultT = None,
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
            + (' [reusable]' if self.reusable else '')
            + (' [skipped]' if self.skipped else '')
            + (f' [error: {self.error_msg}]' if self.error_msg else '')
            + f' meta={self.meta}'
            + ')'
        )
        return res

    def _get(self, key=None) -> str | Path:
        if self.data is None:
            raise ValueError(f'{self.stage}: output data is not available')

        if key is not None:
            if not isinstance(self.data, dict):
                raise ValueError(f'{self.stage}: {self.data} is not a dictionary, can\'t get "{key}"')
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
        if not output.target.get_sequencing_groups():
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
                f'@stage(required_stages=[{stage.__name__}])',
            )

        if stage.__name__ not in self._outputs_by_target_by_stage:
            raise WorkflowError(
                f'No inputs from {stage.__name__} for {self.stage.name} found '
                + 'after skipping targets with missing inputs. '
                + (
                    'Check the logs if all sequencing groups were missing inputs from previous '
                    'stages, and consider changing `workflow/first_stage`'
                    if get_config()['workflow'].get('skip_sgs_with_missing_input')
                    else ''
                ),
            )

        return {trg: fun(result) for trg, result in self._outputs_by_target_by_stage.get(stage.__name__, {}).items()}

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
                f'{self.stage.name}. Is {stage.__name__} in the `required_stages`'
                f'decorator? Available: {self._outputs_by_target_by_stage}',
            )
        if not self._outputs_by_target_by_stage[stage.__name__].get(target.target_id):
            raise StageInputNotFoundError(
                f'Not found output for {target} from stage {stage.__name__}, required for stage {self.stage.name}',
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
        target_sequencing_groups = target.get_sequencing_group_ids()
        for stage_, outputs_by_target in self._outputs_by_target_by_stage.items():
            for target_, output in outputs_by_target.items():
                if output:
                    output_sequencing_groups = output.target.get_sequencing_group_ids()
                    sequencing_groups_intersect = set(target_sequencing_groups) & set(output_sequencing_groups)
                    if sequencing_groups_intersect:
                        for j in output.jobs:
                            assert j, f'Stage: {stage_}, target: {target_}, output: {output}'
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
    i.e. SequencingGroupStage(Stage[SequencingGroup]) should only be able to work on SequencingGroup(Target).
    """

    def __init__(
        self,
        name: str,
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        analysis_type: str | None = None,
        analysis_keys: list[str] | None = None,
        update_analysis_meta: Callable[[str], dict] | None = None,
        tolerate_missing_output: bool = False,
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
        # If `analysis_keys` are defined, it will be used to extract the value for
        # `Analysis.output` if the Stage.expected_outputs() returns a dict.
        self.analysis_keys = analysis_keys
        # if `update_analysis_meta` is defined, it is called on the `Analysis.output`
        # field, and result is merged into the `Analysis.meta` dictionary.
        self.update_analysis_meta = update_analysis_meta

        self.tolerate_missing_output = tolerate_missing_output

        # Populated with the return value of `add_to_the_workflow()`
        self.output_by_target: dict[str, StageOutput | None] = dict()

        self.skipped = skipped
        self.forced = forced or self.name in get_config()['workflow'].get('force_stages', [])
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

    @property
    def analysis_prefix(self) -> Path:
        return get_workflow().analysis_prefix / self.name

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

    def deprecated_queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Queues jobs for each corresponding target, defined by Stage subclass.

        Returns a dictionary of `StageOutput` objects indexed by target unique_id.
        """
        return {}

    @abstractmethod
    def queue_for_multicohort(self, multicohort: MultiCohort) -> dict[str, StageOutput | None]:
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
        data: ExpectedResultT = None,
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
                    assert input_job, f'Input dependency job for stage: {self}, target: {target}'
                    output_job.depends_on(input_job)

        if outputs.error_msg:
            return outputs

        # Adding status reporter jobs
        if self.analysis_type and self.status_reporter and action == Action.QUEUE and outputs.data:
            analysis_outputs: list[str | Path] = []
            if isinstance(outputs.data, dict):
                if not self.analysis_keys:
                    raise WorkflowError(
                        f'Cannot create Analysis: `analysis_keys` '
                        f'must be set with the @stage decorator to select value from '
                        f'the expected_outputs dict: {outputs.data}',
                    )
                if not all(key in outputs.data for key in self.analysis_keys):
                    raise WorkflowError(
                        f'Cannot create Analysis for stage {self.name}: `analysis_keys` '
                        f'"{self.analysis_keys}" is not a subset of the expected_outputs '
                        f'keys {outputs.data.keys()}',
                    )

                for analysis_key in self.analysis_keys:
                    analysis_outputs.append(outputs.data[analysis_key])

            else:
                analysis_outputs.append(outputs.data)

            project_name = None
            if isinstance(target, SequencingGroup):
                project_name = target.dataset.name
            elif isinstance(target, Dataset):
                project_name = target.name
            elif isinstance(target, (Cohort, MultiCohort)):
                project_name = target.analysis_dataset.name

            assert isinstance(project_name, str)

            # bump name to include `-test`
            if get_config()['workflow']['access_level'] == 'test' and 'test' not in project_name:
                project_name = f'{project_name}-test'

            for analysis_output in analysis_outputs:
                if not outputs.jobs:
                    continue

                assert isinstance(analysis_output, (str, Path)), f'{analysis_output} should be a str or Path object'
                if outputs.meta is None:
                    outputs.meta = {}

                self.status_reporter.create_analysis(
                    b=get_batch(),
                    output=str(analysis_output),
                    analysis_type=self.analysis_type,
                    target=target,
                    jobs=outputs.jobs,
                    job_attr=self.get_job_attrs(target) | {'stage': self.name, 'tool': 'metamist'},
                    meta=outputs.meta,
                    update_analysis_meta=self.update_analysis_meta,
                    tolerate_missing_output=self.tolerate_missing_output,
                    project_name=project_name,
                )

        return outputs

    def _get_action(self, target: TargetT) -> Action:
        """
        Based on stage parameters and expected outputs existence, determines what
        to do with the target: queue, skip or reuse, etc...
        """
        if target.forced and not self.skipped:
            logging.info(f'{self.name}: {target} [QUEUE] (target is forced)')
            return Action.QUEUE

        if (d := get_config()['workflow'].get('skip_stages_for_sgs')) and self.name in d:
            skip_targets = d[self.name]
            if target.target_id in skip_targets:
                logging.info(f'{self.name}: {target} [SKIP] (is in workflow/skip_stages_for_sgs)')
                return Action.SKIP

        expected_out = self.expected_outputs(target)
        reusable, first_missing_path = self._is_reusable(expected_out)

        if self.skipped:
            if reusable and not first_missing_path:
                logging.info(f'{self.name}: {target} [REUSE] (stage skipped, and outputs exist)')
                return Action.REUSE
            if get_config()['workflow'].get('skip_sgs_with_missing_input'):
                logging.warning(
                    f'{self.name}: {target} [SKIP] (stage is required, '
                    f'but is marked as "skipped", '
                    f'workflow/skip_sgs_with_missing_input=true '
                    f'and some expected outputs for the target do not exist: '
                    f'{first_missing_path}',
                )
                # `workflow/skip_sgs_with_missing_input` means that we can ignore
                # sgs/datasets that have missing results from skipped stages.
                # This is our case, so indicating that this sg/dataset should
                # be ignored:
                target.active = False
                return Action.SKIP
            if self.name in get_config()['workflow'].get('allow_missing_outputs_for_stages', []):
                logging.info(
                    f'{self.name}: {target} [REUSE] (stage is skipped, some outputs are'
                    f'missing, but stage is listed in '
                    f'workflow/allow_missing_outputs_for_stages)',
                )
                return Action.REUSE
            else:
                raise WorkflowError(
                    f'{self.name}: stage is required, but is skipped, and '
                    f'the following expected outputs for target {target} do not exist: '
                    f'{first_missing_path}',
                )

        if reusable and not first_missing_path:
            if target.forced:
                logging.info(
                    f'{self.name}: {target} [QUEUE] (can reuse, but forcing the target to rerun this stage)',
                )
                return Action.QUEUE
            elif self.forced:
                logging.info(f'{self.name}: {target} [QUEUE] (can reuse, but forcing the stage to rerun)')
                return Action.QUEUE
            else:
                logging.info(f'{self.name}: {target} [REUSE] (expected outputs exist: {expected_out})')
                return Action.REUSE

        logging.info(f'{self.name}: {target} [QUEUE]')
        return Action.QUEUE

    def _is_reusable(self, expected_out: ExpectedResultT) -> tuple[bool, Path | None]:
        """
        Checks if the outputs of prior stages already exist, and can be reused
        Args:
            expected_out (ExpectedResultT): expected outputs of a stage

        Returns:
            tuple[bool, Path | None]:
                bool: True if the outputs can be reused, False otherwise
                Path | None: first missing path, if any
        """
        if self.assume_outputs_exist:
            logging.info(f'Assuming outputs exist. Expected output is {expected_out}')
            return True, None

        if not expected_out:
            # Marking is reusable. If the stage does not naturally produce any outputs,
            # it would still need to create some flag file.
            logging.info('No expected outputs, assuming outputs exist')
            return True, None

        if get_config()['workflow'].get('check_expected_outputs'):
            paths = path_walk(expected_out)
            logging.info(f'Checking if {paths} from expected output {expected_out} exist')
            if not paths:
                logging.info(f'{expected_out} is not reusable. No paths found.')
                return False, None

            if first_missing_path := next((p for p in paths if not exists(p)), None):
                logging.info(f'{expected_out} is not reusable, {first_missing_path} is missing')
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
    analysis_keys: list[str | Path] | None = None,
    update_analysis_meta: Callable[[str], dict] | None = None,
    tolerate_missing_output: bool = False,
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
    class GenotypeSample(SequencingGroupStage):
        def expected_outputs(self, sequencing_group: SequencingGroup):
            ...
        def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
            ...

    @analysis_type: if defined, will be used to create/update `Analysis` entries
        using the status reporter.
    @analysis_keys: is defined, will be used to extract the value for `Analysis.output`
        if the Stage.expected_outputs() returns a dict.
    @update_analysis_meta: if defined, this function is called on the `Analysis.output`
        field, and returns a dictionary to be merged into the `Analysis.meta`
    @tolerate_missing_output: if True, when registering the output of this stage,
        allow for the output file to be missing (only relevant for metamist entry)
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
                analysis_keys=analysis_keys,
                update_analysis_meta=update_analysis_meta,
                skipped=skipped,
                assume_outputs_exist=assume_outputs_exist,
                forced=forced,
                tolerate_missing_output=tolerate_missing_output,
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
    reason: str | None = None,
    assume_outputs_exist: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Decorator on top of `@stage` that sets the `self.skipped` field to True.
    By default, expected outputs of a skipped stage will be checked,
    unless `assume_outputs_exist` is True.

    @skip
    @stage
    class MyStage1(SequencingGroupStage):
        ...

    @skip
    @stage(assume_outputs_exist=True)
    class MyStage2(SequencingGroupStage):
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
    Encapsulates a Hail Batch object, stages, and a cohort of datasets of sequencing groups.
    Responsible for orchestrating stages.
    """

    def __init__(
        self,
        stages: list[StageDecorator] | None = None,
        dry_run: bool | None = None,
    ):
        if _workflow is not None:
            raise ValueError('Workflow already initialised. Use get_workflow() to get the instance')

        self.dry_run = dry_run or get_config(True)['workflow'].get('dry_run')

        analysis_dataset = get_config(True)['workflow']['dataset']
        name = get_config()['workflow'].get('name', analysis_dataset)
        description = get_config()['workflow'].get('description', name)
        self.name = slugify(name)

        self._output_version: str | None = None
        if output_version := get_config()['workflow'].get('output_version'):
            self._output_version = slugify(output_version)

        self.run_timestamp: str = get_config()['workflow'].get('run_timestamp') or timestamp()

        # Description
        if self._output_version:
            description += f': output_version={self._output_version}'
        description += f': run_timestamp={self.run_timestamp}'
        if sequencing_type := get_config()['workflow'].get('sequencing_type'):
            description += f' [{sequencing_type}]'
        if not self.dry_run:
            if ds_set := set(d.name for d in get_multicohort().get_datasets()):
                description += ' ' + ', '.join(sorted(ds_set))
            reset_batch()
            get_batch().name = description

        self.status_reporter = None
        if get_config()['workflow'].get('status_reporter') == 'metamist':
            self.status_reporter = MetamistStatusReporter()
        self._stages: list[StageDecorator] | None = stages
        self.queued_stages: list[Stage] = []

    @property
    def output_version(self) -> str:
        return self._output_version or get_multicohort().alignment_inputs_hash()

    @property
    def analysis_prefix(self) -> Path:
        return self._prefix(category='analysis')

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
        return get_multicohort().analysis_dataset.prefix(category=category) / self.name / self.output_version

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
                        f'{", ".join(stage_names)}',
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
                    logging.info(f'Skipping stage {descendant} (precedes {fs} listed in first_stages)')
                    stages_d[descendant].skipped = True
                for grand_descendant in nx.descendants(graph, descendant):
                    if not stages_d[grand_descendant].assume_outputs_exist:
                        logging.info(
                            f'Not checking expected outputs of not immediately '
                            f'required stage {grand_descendant} (< {descendant} < {fs})',
                        )
                        stages_d[grand_descendant].assume_outputs_exist = True

            for ancestor in nx.ancestors(graph, fs):
                first_stages_keeps.append(ancestor)

        for ls in last_stages:
            # ancestors of this last_stage
            ancestors = nx.ancestors(graph, ls)
            if any(anc in last_stages for anc in ancestors):
                # a downstream stage is also in last_stages, so this is not yet
                # a "real" last stage that we want to run
                continue
            for ancestor in ancestors:
                if stages_d[ancestor].skipped:
                    continue  # already skipped
                logging.info(f'Skipping stage {ancestor} (after last {ls})')
                stages_d[ancestor].skipped = True
                stages_d[ancestor].assume_outputs_exist = True

            for ancestor in nx.descendants(graph, ls):
                last_stages_keeps.append(ancestor)

        for _stage in stages:
            if _stage.name not in last_stages_keeps + first_stages_keeps:
                _stage.skipped = True
                _stage.assume_outputs_exist = True

    @staticmethod
    def _process_only_stages(stages: list[Stage], graph: nx.DiGraph, only_stages: list[str]):
        if not only_stages:
            return

        stages_d = {s.name: s for s in stages}
        stage_names = list(stg.name for stg in stages)
        lower_names = {s.lower() for s in stage_names}

        for s_name in only_stages:
            if s_name.lower() not in lower_names:
                raise WorkflowError(
                    f'Value in workflow/only_stages "{s_name}" must be a stage '
                    f'name or a subset of stages from the available list: '
                    f'{", ".join(stage_names)}',
                )

        # We want to run stages only appearing in only_stages, and check outputs of
        # imediate predecessor stages, but skip everything else.
        required_stages: set[str] = set()
        for os in only_stages:
            rs = nx.descendants_at_distance(graph, os, 1)
            required_stages |= set(rs)

        for stage in stages:
            # Skip stage not in only_stages, and assume outputs exist...
            if stage.name not in only_stages:
                stage.skipped = True
                stage.assume_outputs_exist = True

        # ...unless stage is directly required by any stage in only_stages
        for stage_name in required_stages:
            stages_d[stage_name].assume_outputs_exist = False

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

        # Only allow one of only_stages or first_stages/last_stages as they seem
        # to be mutually exclusive.
        if only_stages and (first_stages or last_stages or skip_stages):
            raise WorkflowError(
                "Workflow config parameter 'only_stages' is incompatible with "
                + "'first_stages', 'last_stages' and/or 'skip_stages'",
            )

        logging.info(f'End stages for the workflow "{self.name}": {[cls.__name__ for cls in requested_stages]}')
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
                logging.info(f'Additional implicit stages: {list(newly_implicitly_added_d.keys())}')
                _stages_d |= newly_implicitly_added_d
            else:
                # No new implicit stages added, so can stop the depth-search here
                break

        # Round 3: set "stage.required_stages" fields to each stage.
        for stg in _stages_d.values():
            stg.required_stages = [
                _stages_d[cls.__name__] for cls in stg.required_stages_classes if cls.__name__ in _stages_d
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
        if first_stages or last_stages:
            logging.info('Applying workflow/first_stages and workflow/last_stages')
            self._process_first_last_stages(stages, dag, first_stages, last_stages)
        elif only_stages:
            logging.info('Applying workflow/only_stages')
            self._process_only_stages(stages, dag, only_stages)

        if not (final_set_of_stages := [s.name for s in stages if not s.skipped]):
            raise WorkflowError('No stages to run')

        logging.info(f'Final list of stages after applying stage configuration options:\n{final_set_of_stages}')

        required_skipped_stages = [s for s in stages if s.skipped]
        if required_skipped_stages:
            logging.info(f'Skipped stages: {", ".join(s.name for s in required_skipped_stages)}')

        # Round 6: actually adding jobs from the stages.
        if not self.dry_run:
            inputs = get_multicohort()  # Would communicate with metamist.
            for i, stg in enumerate(stages):
                logging.info('*' * 60)
                logging.info(f'Stage #{i + 1}: {stg}')
                # NOTE: When manual dataset entry is deprecated, we can remove this
                if isinstance(inputs, Cohort):
                    stg.output_by_target = stg.deprecated_queue_for_cohort(inputs)
                elif isinstance(inputs, MultiCohort):
                    stg.output_by_target = stg.queue_for_multicohort(inputs)
                else:
                    raise WorkflowError(f'Unsupported input type: {inputs}')
                if errors := self._process_stage_errors(stg.output_by_target):
                    raise WorkflowError(f'Stage {stg} failed to queue jobs with errors: ' + '\n'.join(errors))

                logging.info('')

        else:
            self.queued_stages = [stg for stg in _stages_d.values() if not stg.skipped]
            logging.info(f'Queued stages: {self.queued_stages}')

    @staticmethod
    def _process_stage_errors(output_by_target: dict[str, StageOutput | None]) -> list[str]:
        targets_by_error = defaultdict(list)
        for target, output in output_by_target.items():
            if output and output.error_msg:
                targets_by_error[output.error_msg].append(target)
        return [f'{error}: {", ".join(target_ids)}' for error, target_ids in targets_by_error.items()]


class SequencingGroupStage(Stage[SequencingGroup], ABC):
    """
    Sequencing Group level stage.
    """

    @abstractmethod
    def expected_outputs(self, sequencing_group: SequencingGroup) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def deprecated_queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()

        if not (datasets := cohort.get_datasets()):
            logging.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.input_datasets` or `workflow.input_cohorts` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`',
            )
            return output_by_target
        if not cohort.get_sequencing_groups():
            logging.warning(
                f'{len(cohort.get_sequencing_groups())}/'
                f'{len(cohort.get_sequencing_groups(only_active=False))} '
                f'usable (active=True) sequencing groups found. Check logs above for '
                f'possible reasons sequencing groups were skipped (e.g. all sequencing groups ignored '
                f'via `workflow.skip_sgs` in config, or they all missing stage '
                f'inputs and `workflow.skip_sgs_with_missing_input=true` is set)',
            )
            return output_by_target

        for dataset in datasets:
            if not dataset.get_sequencing_groups():
                logging.warning(
                    f'{dataset}: '
                    f'{len(dataset.get_sequencing_groups())}/'
                    f'{len(dataset.get_sequencing_groups(only_active=False))} '
                    f'usable (active=True) sequencing groups found. Check logs above for '
                    f'possible reasons sequencing groups were skipped (e.g. all sequencing groups ignored '
                    f'via `workflow.skip_sgs` in config, or they all missing stage '
                    f'inputs and `workflow.skip_sgs_with_missing_input=true` is set)',
                )
                continue

            logging.info(f'Dataset {dataset}:')
            # collect all expected outputs across all samples
            # find all directories which will be checked
            # list outputs in advance
            all_outputs: set[Path] = set()
            for sequencing_group in dataset.get_sequencing_groups():
                all_outputs = path_walk(self.expected_outputs(sequencing_group), all_outputs)

            # evaluate_stuff en masse
            for sequencing_group in dataset.get_sequencing_groups():
                action = self._get_action(sequencing_group)
                output_by_target[sequencing_group.target_id] = self._queue_jobs_with_checks(sequencing_group, action)

        return output_by_target

    def queue_for_multicohort(self, multicohort: MultiCohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()

        if not (cohorts := multicohort.get_cohorts()):
            logging.warning(
                f'{len(multicohort.get_cohorts())}/'
                f'{len(multicohort.get_cohorts(only_active=False))} '
                f'usable (active=True) cohorts found in the multicohort. Check that '
                f'`workflow.input_cohorts` is provided, and not all cohorts are skipped ',
            )
            return output_by_target

        for cohort in cohorts:
            if not (datasets := cohort.get_datasets()):
                logging.warning(
                    f'{len(cohort.get_datasets())}/'
                    f'{len(cohort.get_datasets(only_active=False))} '
                    f'usable (active=True) datasets found in the cohort. Check that '
                    f'`workflow.input_datasets` or `workflow.input_cohorts` is provided, and each cohort is populated in metamist.',
                )
                continue

            for dataset in datasets:
                if not dataset.get_sequencing_groups():
                    logging.warning(
                        f'{dataset}: '
                        f'{len(dataset.get_sequencing_groups())}/'
                        f'{len(dataset.get_sequencing_groups(only_active=False))} '
                        f'usable (active=True) sequencing groups found. Check logs above for '
                        f'possible reasons sequencing groups were skipped (e.g. all sequencing groups ignored '
                        f'via `workflow.skip_sgs` in config, or they all missing stage '
                        f'inputs and `workflow.skip_sgs_with_missing_input=true` is set)',
                    )
                    continue

                logging.info(f'Dataset {dataset}:')
                # collect all expected outputs across all samples
                # find all directories which will be checked
                # list outputs in advance
                all_outputs: set[Path] = set()
                for sequencing_group in dataset.get_sequencing_groups():
                    all_outputs = path_walk(self.expected_outputs(sequencing_group), all_outputs)

                # evaluate_stuff en masse
                for sequencing_group in dataset.get_sequencing_groups():
                    action = self._get_action(sequencing_group)
                    output_by_target[sequencing_group.target_id] = self._queue_jobs_with_checks(
                        sequencing_group,
                        action,
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

    def deprecated_queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        if not (datasets := cohort.get_datasets()):
            logging.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.input_datasets` or `workflow.input_cohorts` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`',
            )
            return output_by_target
        for dataset_i, dataset in enumerate(datasets):
            action = self._get_action(dataset)
            logging.info(f'{self.name}: #{dataset_i + 1}/{dataset} [{action.name}]')
            output_by_target[dataset.target_id] = self._queue_jobs_with_checks(dataset, action)
        return output_by_target

    def queue_for_multicohort(self, multicohort: MultiCohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        for cohort in multicohort.get_cohorts():
            for dataset_i, dataset in enumerate(cohort.get_datasets()):
                action = self._get_action(dataset)
                logging.info(f'{self.name}: #{dataset_i + 1}/{dataset} [{action.name}]')
                output_by_target[dataset.target_id] = self._queue_jobs_with_checks(dataset, action)
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

    def deprecated_queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        return {cohort.target_id: self._queue_jobs_with_checks(cohort)}

    def queue_for_multicohort(self, multicohort: MultiCohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        for cohort in multicohort.get_cohorts():
            action = self._get_action(cohort)
            logging.info(f'{self.name}: {cohort} [{action.name}]')
            output_by_target[cohort.target_id] = self._queue_jobs_with_checks(cohort, action)
        return output_by_target


class MultiCohortStage(Stage, ABC):
    """
    MultiCohort-level stage (all datasets of a workflow run).
    """

    @abstractmethod
    def expected_outputs(self, multicohort: MultiCohort) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """
        pass

    @abstractmethod
    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_multicohort(self, multicohort: MultiCohort) -> dict[str, StageOutput | None]:
        """
        Plug the stage into the workflow.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        action = self._get_action(multicohort)
        logging.info(f'{self.name}: {multicohort} [{action.name}]')
        output_by_target[multicohort.target_id] = self._queue_jobs_with_checks(multicohort, action)
        return output_by_target
