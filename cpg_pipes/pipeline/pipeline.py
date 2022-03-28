"""
Provides `Pipeline` class and a `@stage` decorator that allow to define pipeline
stages and plug them together.

A Stage describes the job that are added to Hail Batch, and outputs that are expected
to be produced. Each stage acts on a `Target`, which can be of a different level:
a `Sample`, a `Dataset`, a `Cohort` (= all input datasets combined). Pipeline plugs
stages together by resolving dependencies between different levels accordingly.

Examples of pipelines can be found in the `pipelines/` folder in the repository root.
"""
import functools
import logging
import pathlib
import shutil
import tempfile
import time
from abc import ABC, abstractmethod
from typing import Callable, cast, Union, TypeVar, Generic, Any, Optional, Type

from cloudpathlib import CloudPath
import hailtop.batch as hb
from hailtop.batch.job import Job

from .exceptions import PipelineError
from .targets import Target, Dataset, Sample, Cohort
from .. import Path, to_path
from ..providers import (
    Cloud,
    Namespace,
    StoragePolicy,
    StatusReporterType,
    InputProviderType,
)
from ..providers.inputs import InputProvider, CsvInputProvider
from ..providers.status import StatusReporter
from ..providers.cpg import (
    CpgStorageProvider,
    SmdbStatusReporter,
    SmdbInputProvider,
    SMDB,
)
from ..refdata import RefData
from ..utils import exists
from ..hb.batch import Batch, setup_batch, get_billing_project
from ..hb.prev_job import PrevJob

logger = logging.getLogger(__file__)


# list for the case when the user doesn't pass them explicitly with set_stages()
_ALL_DEFINED_STAGES = []


StageDecorator = Callable[..., 'Stage']

# Type variable to make sure a Stage subclass always matches the
# correspondinng Target subclass
TargetT = TypeVar('TargetT', bound=Target)

ExpectedResultT = Union[Path, dict[str, Path], None]

StageOutputData = Union[
    Path, hb.Resource, dict[str, Path], dict[str, hb.Resource]
]


# noinspection PyShadowingNames
class StageOutput:
    """
    Represents a result of a specific stage, which was run on a specific target.
    Can be a file path, or a Hail Batch Resource. Optionally wrapped in a dict.
    """
    def __init__(
        self,
        data: StageOutputData,
        stage: 'Stage',
        target: 'Target',
        jobs: list[Job] | None = None,
    ):
        self.data = data
        self.stage = stage
        self.target = target
        self.jobs: list[Job] = jobs or []

    def __repr__(self) -> str:
        return f'result {self.data} for target {self.target}, stage {self.stage}'

    def as_path_or_resource(self, id=None) -> Path | hb.Resource:
        """
        Cast the result to Union[str, hb.Resource], error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        if id is not None:
            if not isinstance(self.data, dict):
                raise ValueError(
                    f'{self.data} is not a dictionary, can\'t get "{id}".'
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
    Represents an input for a stage run. It wraps the outputs of all required upstream stages
    for corresponding targets (e.g. all GVCFs from a HaploytypeCallerStage
    for a JointCallingStage, along with Hail Batch jobs).

    An object of this class is passed to the public `queue_jobs` method of a Stage, 
    and can be used to query dependency files and jobs.
    """

    def __init__(self, stage: 'Stage'):
        self.stage = stage
        self._results_by_target_by_stage: dict[str, dict[str, StageOutput]] = {}
        self._jobs: list[Job] = []

    def add_other_stage_output(self, output: StageOutput):
        """
        Add output from another stage run
        """
        if output.target.active:
            stage_name = output.stage.name
            target_id = output.target.target_id

            if stage_name not in self._results_by_target_by_stage:
                self._results_by_target_by_stage[stage_name] = dict()
            self._results_by_target_by_stage[stage_name][target_id] = output

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

        if stage.__name__ not in self._results_by_target_by_stage:
            raise PipelineError(
                f'No inputs from {stage.__name__} for {self.stage.name} found '
                'after skipping targets with missing inputs. ' +
                ('Check the logs if all samples were missing inputs from previous '
                 'stages, and consider changing --first-stage'
                 if self.stage.skip_samples_with_missing_input else '')
            )

        return {
            trg: fun(result)
            for trg, result
            in self._results_by_target_by_stage.get(stage.__name__, {}).items()
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

    def as_dict_by_target(
        self, stage: StageDecorator
    ) -> dict[str, dict[str, Path]]:
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
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
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
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_resource(id)

    def as_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dict of files or Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_dict()

    def as_path_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, Path]:
        """
        Get a dict of files for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_path_dict()

    def as_resource_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, hb.Resource]:
        """
        Get a dict of  Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_resource_dict()

    def get_jobs(self) -> list[Job]:
        """
        Build a list of hail batch dependencies from all stages and targets
        """
        all_jobs = []
        for _, results_by_target in self._results_by_target_by_stage.items():
            for _, results in results_by_target.items():
                all_jobs.extend(results.jobs)
        return all_jobs


class Stage(Generic[TargetT], ABC):
    """
    Abstract class for a pipeline stage.
    """

    def __init__(
        self,
        name: str,
        batch: Batch,
        refs: RefData,
        hail_billing_project: str,
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        analysis_type: str | None = None,
        skipped: bool = False,
        required: bool | None = None,
        assume_outputs_exist: bool = False,
        forced: bool = False,
        skip_samples_with_missing_input: bool = False,
        check_expected_outputs: bool = False,
        check_intermediates: bool = False,
        status_reporter: StatusReporter | None = None,
        pipeline_config: dict[str, Any] | None = None,
        hail_bucket: Path | None = None,
    ):
        """
        @param name: name of the stage
        @param batch: Hail Batch object
        @param refs: reference data for bioinformatics
        @param hail_billing_project: Hail Batch billing project
        @param required_stages: list of stage classes that this stage requires
        @param analysis_type: if defined, will query the SMDB Analysis entries 
            of this type
        @param skipped: means that the stage is skipped and self.queue_jobs()
            won't run. The other stages if depend on it can aassume that that 
            self.expected_outputs() returns existing files and target.ouptut_by_stage 
            will be populated.
        @param required: means that the self.expected_output() results are 
            required for another active stage, even if the stage was skipped.
        @param assume_outputs_exist: for skipped but required stages, 
            the self.expected_outputs() output will still be checked for existence. 
            This option makes the downstream stages assume that the output exist.
        @param forced: run self.queue_jobs(), even if we can reuse the 
            self.expected_output().
        @param hail_bucket: Hail Batch bucket
        """
        self._name = name
        self.b = batch
        self.refs = refs

        self.skip_samples_with_missing_input = skip_samples_with_missing_input
        self.check_expected_outputs = check_expected_outputs
        self.check_intermediates = check_intermediates

        self.required_stages_classes: list[StageDecorator] = []
        if required_stages:
            if isinstance(required_stages, list):
                self.required_stages_classes.extend(required_stages)
            else:
                self.required_stages_classes.append(required_stages)

        # Populated in pipeline.run(), after we know all stages
        self.required_stages: list[Stage] = []

        self.status_reporter = status_reporter
        # If analysis type is defined, it will be used to update analysis status,
        # as well as find and reuse existing outputs from the status reporter
        self.analysis_type = analysis_type

        # Populated with the return value of `add_to_the_pipeline()`
        self.output_by_target: dict[str, StageOutput] = dict()

        self.skipped = skipped
        self.required = required if required is not None else not skipped
        self.forced = forced
        self.assume_outputs_exist = assume_outputs_exist

        self.pipeline_config = pipeline_config or {}
        
        self.hail_billing_project = hail_billing_project
        self.hail_bucket = hail_bucket

    @property
    def name(self):
        """
        Stage name (unique and descriptive stage)
        """
        return self._name

    @abstractmethod
    def queue_jobs(self, target: TargetT, inputs: StageInput) -> StageOutput:
        """
        Implements logic of the Stage: creates Batch jobs that do the processing.
        Assumes that all the household work is done: checking missing inputs
        from requried stages, checking for possible reuse of existing outputs.
        """

    @abstractmethod
    def expected_outputs(self, target: TargetT) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `queue_jobs()`.

        Can be a str or a AnyPath object, or a dictionary of str/Path objects.
        """

    @abstractmethod
    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> dict[str, StageOutput]:
        """
        Calls `output = pipeline.add_for_target(target)` on each target, 
        which itself calls `output = queue_jobs(target, input)`, making sure to
        construct the correct `input`.

        Returns a dictionary of StageOutput, indexed by target unique_id.
        """

    def make_outputs(
        self,
        target: TargetT,
        data: StageOutputData | str | dict[str, str] | None = None,
        jobs: list[Job] | Job | None = None
    ) -> StageOutput:
        """
        Builds a StageDeps object to return from a stage's queue_jobs()
        """
        # Converting str into Path objects.
        path_data: StageOutputData
        if isinstance(data, dict):
            path_data = {k: to_path(v) for k, v in data.items()}
        elif data is not None:
            path_data = to_path(data)
        else:
            path_data = data
        jobs = [jobs] if isinstance(jobs, Job) else jobs
        # Adding status reporter jobs
        if self.analysis_type:
            if (
                isinstance(path_data, hb.Resource) or (
                    isinstance(path_data, dict) and 
                    any(isinstance(d, hb.Resource) for k, d in path_data.items())
                )
            ):
                raise PipelineError(
                    'Cannot use hb.Resource objects with status reporter. '
                    'Only supported Path objects and dicts of Path objects'
                )
            if self.status_reporter:
                self.status_reporter.add_updaters_jobs(
                    b=self.b,
                    output=path_data,
                    analysis_type=self.analysis_type,
                    target=target,
                    jobs=jobs,
                )
        return StageOutput(stage=self, target=target, data=data, jobs=jobs)

    def _make_inputs(self) -> StageInput:
        """
        Collects outputs from all dependencies and create input for this stage
        """
        inputs = StageInput(self)
        for prev_stage in self.required_stages:
            for _, stage_output in prev_stage.output_by_target.items():
                inputs.add_other_stage_output(stage_output)
        return inputs

    def _queue_jobs_with_checks(self, target: TargetT) -> StageOutput:
        """
        Constructs `inputs` and calls the public `output = queue_jobs(target, input)`.
        Performs the checks like possibility to reuse existing jobs results,
        or if required dependencies are missing.
        """
        if not self.skipped:
            inputs = self._make_inputs()

            reusable_paths = self._try_get_reusable_paths(target)
            if reusable_paths:
                if target.forced:
                    logger.info(
                        f'{self.name}: can reuse, but forcing the target '
                        f'{target} to rerun this stage'
                    )
                    outputs = self.queue_jobs(target, inputs)
                elif self.forced:
                    logger.info(
                        f'{self.name}: can reuse, but forcing the stage '
                        f'to rerun, target={target}'
                    )
                    outputs = self.queue_jobs(target, inputs)
                else:
                    logger.info(f'{self.name}: reusing results for {target}')
                    outputs = self._queue_reuse_job(target, reusable_paths)
            else:
                logger.info(f'{self.name}: adding jobs for {target}')
                outputs = self.queue_jobs(target, inputs)

            for j in outputs.jobs:
                j.depends_on(*inputs.get_jobs())

            return outputs

        elif self.required:
            reusable_paths = self._try_get_reusable_paths(target)
            if not reusable_paths:
                raise ValueError(
                    f'Stage {self.name} is required, but is skipped, and '
                    f'expected outputs for target {target} do not exist'
                )
            else:
                return self.make_outputs(target=target, data=reusable_paths)

        else:
            # Stage is not needed, returning empty outputs
            return self.make_outputs(target=target)

    def _try_get_reusable_paths(self, target: TargetT) -> ExpectedResultT:
        """
        Returns outputs that can be reused for the stage for the target,
        or None of none can be reused
        """
        expected_output = self.expected_outputs(target)

        if not expected_output:
            return None
        elif not self.check_expected_outputs:
            if self.assume_outputs_exist:
                # Do not check the files' existence, trust they exist:
                return expected_output
            else:
                # Do not check the files' existence, assume they don't exist:
                return None
        elif not self.required:
            # This stage is not required, so can just assume outputs exist:
            return expected_output
        else:
            # Checking that expected output exists:
            paths: list[Path]
            if isinstance(expected_output, dict):
                paths = [v for k, v in expected_output.items()]
            else:
                paths = [expected_output]
            if not all(exists(path) for path in paths):
                return None
            return expected_output

    def _queue_reuse_job(
        self,
        target: TargetT,
        found_paths: Path | dict[str, Path]
    ) -> StageOutput:
        """
        Queues a [reuse] Job
        """
        return self.make_outputs(
            target=target,
            data=found_paths,
            jobs=[self.b.new_job(f'{self.name} [reuse]', target.get_job_attrs())]
        )


def stage(
    cls: Optional[Type['Stage']] = None, 
    *,
    analysis_type: str | None = None, 
    required_stages: list[StageDecorator] | StageDecorator | None = None,
    skipped: bool = False,
    required: bool | None = None,
    assume_outputs_exist: bool = False,
    forced: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with an optional argument.
    The goal is to allow cleaner defining of custom pipeline stages, without
    requiring to implement constructor. E.g.

    @stage(analysis_type='gvcf', required_stages=CramStage)
    class GvcfStage(SampleStage):
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
                required_stages=required_stages,
                analysis_type=analysis_type,
                skipped=skipped,
                required=required,
                assume_outputs_exist=assume_outputs_exist, 
                forced=forced,
                skip_samples_with_missing_input=pipeline.skip_samples_with_missing_input,
                check_expected_outputs=pipeline.check_expected_outputs,
                check_intermediates=pipeline.check_intermediates,
                pipeline_config=pipeline.config,
                refs=pipeline.refs,
                hail_billing_project=pipeline.hail_billing_project,
                hail_bucket=pipeline.hail_bucket,
            )
        # We record each initialised Stage subclass, so we know the default stage
        # list for the case when the user doesn't pass them explicitly with set_stages()
        _ALL_DEFINED_STAGES.append(wrapper_stage)
        return wrapper_stage

    if cls is None:
        return decorator_stage
    else:
        return decorator_stage(cls)


def skip(
    _fun: Optional[StageDecorator] = None, 
    *,
    assume_outputs_exist: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Decorator on top of `@stage` that sets the `self.skipped` field to True

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
    Represents a Pipeline, and incapulates a Hail Batch object, stages, 
    and a cohort of datasets of samples.
    """
    def __init__(
        self,
        analysis_dataset: str,
        name: str,
        description: str,
        namespace: Namespace,
        storage_policy: StoragePolicy = StoragePolicy.CPG,
        cloud: Cloud = Cloud.GS,
        status_reporter_type: StatusReporterType = StatusReporterType.NONE,
        input_provider_type: InputProviderType = InputProviderType.NONE,
        input_csv: str | None = None,
        stages_in_order: list[StageDecorator] | None = None,
        keep_scratch: bool = True,
        dry_run: bool = False,
        version: str | None = None,
        previous_batch_tsv_path: Path | None = None,
        previous_batch_id: str | None = None,
        skip_samples_with_missing_input: bool = False,
        check_intermediates: bool = True,
        check_expected_outputs: bool = True,
        first_stage: str | None = None,
        last_stage: str | None = None,
        config: dict | None = None,
        input_datasets: list[str] | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        force_samples: list[str] | None = None,
        ped_files: list[Path] | None = None,
        local_dir: Path | None = None,
    ):
        if stages_in_order and not input_datasets:
            raise ValueError(
                'Datasets must be populated before adding stages. '
                'Provide `input_datasets`, or omit `stages_in_order` and call '
                'pipeline.set_stages(stages_in_order) later.'
            )

        self.name = name
        self.version = version or time.strftime('%Y%m%d-%H%M%S')
        self.namespace = namespace

        if storage_policy == StoragePolicy.CPG:
            storage_provider = CpgStorageProvider(cloud)
        else:
            raise PipelineError(f'Unsupported storage policy {storage_policy}')

        self.cohort = Cohort(
            self.name,
            analysis_dataset_name=analysis_dataset,
            namespace=namespace,
            storage_provider=storage_provider
        )
        self.refs = RefData(storage_provider.get_ref_bucket())

        self.check_intermediates = check_intermediates
        self.check_expected_outputs = check_expected_outputs
        self.skip_samples_with_missing_input = skip_samples_with_missing_input
        self.first_stage = first_stage
        self.last_stage = last_stage

        self.tmp_bucket = to_path(
            self.cohort.analysis_dataset.get_tmp_bucket(
                version=f'{self.name}/{self.version}'
            )
        )
        self.keep_scratch = keep_scratch
        self.dry_run: bool = dry_run

        if local_dir:
            self.local_dir = to_path(local_dir)
            self.local_tmp_dir = to_path(tempfile.mkdtemp(dir=local_dir))
        else:
            self.local_dir = self.local_tmp_dir = to_path(tempfile.mkdtemp())

        self.config = config or {}

        if version:
            description += f' {version}'
        if input_datasets:
            description += ': ' + ', '.join(input_datasets)

        self.hail_billing_project = get_billing_project(
            self.cohort.analysis_dataset.stack
        )
        self.hail_bucket = self.tmp_bucket / 'hail'
        self.b: Batch = setup_batch(
            description=description, 
            billing_project=self.hail_billing_project,
            hail_bucket=self.hail_bucket,
        )
        self.prev_batch_jobs = dict()
        if previous_batch_tsv_path is not None:
            assert previous_batch_id is not None
            self.prev_batch_jobs = PrevJob.parse(
                previous_batch_tsv_path,
                previous_batch_id,
                hail_bucket=self.hail_bucket,
            )

        self.status_reporter: StatusReporter | None = None
        input_provider: InputProvider | None = None
        if (
            input_provider_type == InputProviderType.SMDB 
            or status_reporter_type == StatusReporterType.SMDB
        ):
            smdb = SMDB(self.cohort.analysis_dataset.name)
            if status_reporter_type == StatusReporterType.SMDB:
                self.status_reporter = SmdbStatusReporter(smdb)
            if input_provider_type == InputProviderType.SMDB:
                input_provider = SmdbInputProvider(smdb)

        if input_provider_type == InputProviderType.CSV:
            if not input_csv:
                raise PipelineError(
                    f'input_csv (--input-csv) should be provided '
                    f'with input_provider_type=InputProviderType.CSV '
                    f'(--input-provider {InputProviderType.CSV.value})'
                )
            input_provider = CsvInputProvider(to_path(input_csv).open())

        if input_datasets and input_provider:
            input_provider.populate_cohort(
                cohort=self.cohort,
                dataset_names=input_datasets,
                skip_samples=skip_samples,
                only_samples=only_samples,
                ped_files=ped_files,
            )

        if force_samples:
            for s in self.cohort.get_samples():
                if s.id in force_samples:
                    logger.info(
                        f'Force rerunning sample {s.id} even if its outputs exist'
                    )
                    s.forced = True

        self._stages_dict: dict[str, Stage] = dict()
        self._stages_in_order: list[StageDecorator] = (
            stages_in_order or _ALL_DEFINED_STAGES
        )

    def submit_batch(
        self, 
        dry_run: bool | None = None,
        wait: bool | None = False,
    ) -> Any | None:
        """
        Submits Hail Batch jobs.
        """
        if dry_run is None:
            dry_run = self.dry_run
        
        self.set_stages(self._stages_in_order)

        if self.b:
            logger.info(f'Will submit {self.b.total_job_num} jobs:')
            for label, stat in self.b.labelled_jobs.items():
                logger.info(
                    f'  {label}: {stat["job_n"]} for {len(stat["samples"])} samples'
                )
            logger.info(f'  Other jobs: {self.b.other_job_num}')

            return self.b.run(
                dry_run=dry_run,
                delete_scratch_on_exit=not self.keep_scratch,
                wait=wait,
            )
        shutil.rmtree(self.local_tmp_dir)
        return None

    def _validate_first_last_stage(self) -> tuple[int | None, int | None]:
        # Validating the --first-stage and --last-stage parameters
        stage_names = list(self._stages_dict.keys())
        lower_stage_names = [s.lower() for s in stage_names]
        first_stage_num = None
        if self.first_stage:
            if self.first_stage.lower() not in lower_stage_names:
                logger.critical(
                    f'Value for --first-stage {self.first_stage} '
                    f'not found in available stages: {", ".join(stage_names)}'
                )
            first_stage_num = lower_stage_names.index(self.first_stage.lower())
        last_stage_num = None
        if self.last_stage:
            if self.last_stage.lower() not in lower_stage_names:
                logger.critical(
                    f'Value for --last-stage {self.last_stage} '
                    f'not found in available stages: {", ".join(stage_names)}'
                )
            last_stage_num = lower_stage_names.index(self.last_stage.lower())
        return first_stage_num, last_stage_num

    def set_stages(self, stages_classes: list[StageDecorator]):
        """
        Iterate over stages and call add_to_the_pipeline() on each.
        Effectively creates all Hail Batch jobs through Stage.queue_jobs().
        """
        # Initializing stage objects
        stages = [cls(self) for cls in stages_classes]
        for stage_ in stages:
            if stage_.name in self._stages_dict:
                raise ValueError(
                    f'Stage {stage_.name} is already defined. Check your '
                    f'list for duplicates: {", ".join(s.name for s in stages)}'
                )
            self._stages_dict[stage_.name] = stage_

        first_stage_num, last_stage_num = self._validate_first_last_stage()

        # First round - checking which stages we require, even if they are skipped.
        # If there are required stages that are not defined, we defined them
        # into `additional_stages` as `skipped`.
        additional_stages_dict: dict[str, Stage] = dict()
        for i, (stage_name, stage_) in enumerate(self._stages_dict.items()):
            if first_stage_num is not None and i < first_stage_num:
                stage_.skipped = True
                stage_.required = False
                logger.info(f'Skipping stage {stage_name}')
                continue

            if last_stage_num is not None and i > last_stage_num:
                stage_.skipped = True
                stage_.required = False
                continue

            for reqcls in stage_.required_stages_classes:
                if reqcls.__name__ in self._stages_dict:
                    reqstage = self._stages_dict[reqcls.__name__]
                elif reqcls.__name__ in additional_stages_dict:
                    reqstage = additional_stages_dict[reqcls.__name__]
                else:
                    # stage is not initialized, so we defining it as skipped
                    reqstage = reqcls(self)
                    reqstage.skipped = True
                    additional_stages_dict[reqcls.__name__] = reqstage
                reqstage.required = True
                if reqstage.skipped:
                    logger.info(
                        f'Stage {reqstage.name} is skipped, '
                        f'but the output will be required for the stage {stage_.name}'
                    )
                stage_.required_stages.append(reqstage)
        
        if additional_stages_dict:
            logger.info(
                f'Additional required stages added as "skipped": '
                f'{additional_stages_dict.keys()}'
            )
            for name, stage_ in self._stages_dict.items():
                additional_stages_dict[name] = stage_    
            self._stages_dict = additional_stages_dict

        logger.info(f'Setting stages: {", ".join(s.name for s in stages if not s.skipped)}')
        required_skipped_stages = [s for s in stages if s.skipped and s.required]
        if required_skipped_stages:
            logger.info(
                f'Skipped stages that are used to get inputs: '
                f'{", ".join(s.name for s in required_skipped_stages)}'
            )

        # Second round - actually adding jobs from the stages.
        for i, (_, stage_) in enumerate(self._stages_dict.items()):
            if not stage_.skipped:
                logger.info(f'*' * 60)
                logger.info(f'Stage {stage_.name}')

            if stage_.required:
                logger.info(f'Adding jobs for stage {stage_.name}')
                stage_.output_by_target = stage_.add_to_the_pipeline(self)

            if not stage_.skipped:
                logger.info(f'')
                if last_stage_num and i >= last_stage_num:
                    logger.info(f'Last stage is {stage_.name}, stopping here')
                    break

    def get_datasets(self, only_active: bool = True) -> list[Dataset]:
        """
        Thin wrapper around corresponding Cohort method.
        """
        return self.cohort.get_datasets(only_active=only_active)

    def add_dataset(self, name: str) -> Dataset:
        """
        Thin wrapper around corresponding Cohort method.
        """
        return self.cohort.add_dataset(name)

    def get_all_samples(self, only_active: bool = True) -> list[Sample]:
        """
        Thin wrapper around corresponding Cohort method.
        """
        return self.cohort.get_samples(only_active=only_active)

    def get_all_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Thin wrapper around corresponding Cohort method.
        """
        return self.cohort.get_sample_ids(only_active=only_active)
