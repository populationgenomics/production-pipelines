"""
Stage classes.
"""

import logging
import pathlib
from abc import ABC, abstractmethod
from typing import Callable, cast, Union, TypeVar, Generic, Any

import hailtop.batch as hb
from cloudpathlib import CloudPath

from cpg_pipes.storage import Path, to_path
from hailtop.batch.job import Job

from cpg_pipes.buckets import exists
from cpg_pipes.hb.batch import Batch
from cpg_pipes.pipeline.analysis import AnalysisType
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.exceptions import PipelineError
from cpg_pipes.cpg.smdb import SMDB
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.dataset import Sample, Cohort

logger = logging.getLogger(__file__)


# Type variable to make sure a Stage subclass always matches the
# correspondinng Target subclass
TargetT = TypeVar('TargetT', bound=Target)

ExpectedResultT = Union[Path, dict[str, Path], None]

StageOutputData = Union[
    Path, hb.Resource, dict[str, Path], dict[str, hb.Resource]
]


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


StageDecorator = Callable[..., 'Stage']


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
        Get as a dictoinary of files/resources for a specific stage, indexed by target
        """
        return self._each(fun=(lambda r: r.as_dict()), stage=stage)

    def as_resource_dict_by_target(
        self,
        stage: StageDecorator,
    ) -> dict[str, dict[str, hb.Resource]]:
        """
        Get a dictoinary of resources for a specific stage, and indexed by target
        """
        return self._each(fun=(lambda r: r.as_resource_dict()), stage=stage)

    def as_path_dict_by_target(
        self,
        stage: StageDecorator,
    ) -> dict[str, dict[str, Path]]:
        """
        Get a dictoinary of paths for a specific stage, and indexed by target
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
        Get a dictoinary of files or Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_dict()

    def as_path_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, Path]:
        """
        Get a dictoinary of files for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.target_id]
        return res.as_path_dict()

    def as_resource_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, hb.Resource]:
        """
        Get a dictoinary of  Resources for a specific target and stage
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
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        sm_analysis_type: AnalysisType | None = None,
        skipped: bool = False,
        required: bool | None = None,
        assume_results_exist: bool = False,
        forced: bool = False,
        validate_smdb_analyses: bool = False,
        skip_samples_with_missing_input: bool = False,
        check_expected_outputs: bool = False,
        check_intermediates: bool = False,
        smdb: SMDB | None = None,
        pipeline_config: dict[str, Any] | None = None,
        dry_run: bool = False,
    ):        
        """
        @param required_stages: list of stage classes that this stage requires
        @param sm_analysis_type: if defined, will query the SMDB Analysis entries 
            of this type
        @param skipped: means that the stage is skipped and self.queue_jobs()
            won't run. The other stages if depend on it can aassume that that 
            self.expected_result() returns existing files and target.ouptut_by_stage 
            will be populated.
        @param required: means that the self.expected_output() results are 
            required for another active stage, even if the stage was skipped.
        @param assume_results_exist: for skipped but required stages, 
            the self.expected_result() output will still be checked for existence. 
            This option makes the downstream stages assume that the output exist.
        @param forced: run self.queue_jobs(), even if we can reuse the 
            self.expected_output().
        """
        self._name = name
        self.b = batch
        self.dry_run = dry_run
        
        self.validate_smdb_analyses = validate_smdb_analyses
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

        self.smdb = smdb
        # If analysis type is defined, it will be used to find and reuse existing
        # outputs from the SMDB
        self.analysis_type = sm_analysis_type

        # Populated with the return value of `add_to_the_pipeline()`
        self.output_by_target: dict[str, StageOutput] = dict()

        self.skipped = skipped
        self.required = required if required is not None else not skipped
        self.forced = forced
        self.assume_results_exist = assume_results_exist
        
        self.pipeline_config = pipeline_config or {}

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
        from requried stages, checking the SMDB, checking for possible reuse of 
        existing outputs.
        """

    @abstractmethod
    def expected_result(self, target: TargetT) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `queue_jobs()`.

        Can be a str or a AnyPath object, or a dictionary of str/AnyPath objects.
        """

    @abstractmethod
    def add_to_the_pipeline(self, pipeline) -> dict[str, StageOutput]:
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
        # Converting str into AnyPath objects.
        if isinstance(data, dict):
            data = {k: to_path(v) for k, v in data.items()}
        elif data is not None:
            data = to_path(data)
        jobs = [jobs] if isinstance(jobs, Job) else jobs
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
            reusable_paths = self._try_get_reusable_paths(target)
            if reusable_paths:
                if target.forced:
                    logger.info(
                        f'{self.name}: can reuse, but forcing the target '
                        f'{target.target_id} to rerun this stage'
                    )
                    return self.queue_jobs(target, self._make_inputs())
                elif self.forced:
                    logger.info(
                        f'{self.name}: can reuse, but forcing the stage '
                        f'to rerun, target={target.target_id}'
                    )
                    return self.queue_jobs(target, self._make_inputs())
                else:
                    logger.info(f'{self.name}: reusing results for {target.target_id}')
                    return self._queue_reuse_job(target, reusable_paths)
            else:
                logger.info(f'{self.name}: adding jobs for {target.target_id}')
                return self.queue_jobs(target, self._make_inputs())

        elif self.required:
            reusable_paths = self._try_get_reusable_paths(target)
            if not reusable_paths:
                raise ValueError(
                    f'Stage {self.name} is required, but is skipped, and '
                    f'expected outputs for target {target.target_id} do not exist.)'
                )
            else:
                return self.make_outputs(target=target, data=reusable_paths)

        else:
            # Stage is not needed, returning empty outputs
            return self.make_outputs(target=target)

    def _try_get_reusable_paths(self, target: TargetT) -> dict[str, Path] | Path | None:
        """
        Returns outputs that can be reused for the stage for the target,
        or None of none can be reused
        """
        if self.dry_run:
            return None

        expected_output = self.expected_result(target)

        # Converting all str into AnyPath
        expected_paths: dict[str, Path] | Path | None = None
        if isinstance(expected_output, dict):
            expected_paths = {k: to_path(v) for k, v in expected_output.items()}
        elif expected_output is not None:
            expected_paths = to_path(expected_output)  # type: ignore

        if self.analysis_type is not None:
            # Checking the SMDB
            if not expected_output:
                raise ValueError(
                    f'expected_result() returned None, but must return str '
                    f'for a stage with analysis_type: {self.name} '
                    f'on {target.target_id}, analysis_type={self.analysis_type}'
                )

            if isinstance(expected_paths, dict):
                raise ValueError(
                    f'expected_result() returns a dict, won\'t check the SMDB for '
                    f'{self.name} on {target.target_id}'
                )
            validate = self.validate_smdb_analyses
            if self.skipped and (
                self.assume_results_exist or self.skip_samples_with_missing_input
            ):
                validate = False

            # found_path would be analysis.output from DB if it exists; if it doesn't,
            # and validate_smdb_analyses=True, expected_paths file would be checked:
            found_path = None
            if self.smdb is not None:
                found_path = self._try_get_smdb_analysis(
                    self.smdb, target, cast(Path, expected_paths), validate=validate
                )
            return found_path

        elif not expected_paths:
            return None

        elif (
            not self.check_expected_outputs or
            not self.required or
            self.assume_results_exist
        ):
            # Do not need to check file existence, trust it exists:
            return expected_paths

        else:
            # Checking that expected output exists:
            paths: list[Path]
            if isinstance(expected_paths, dict):
                paths = list(expected_paths.values())
            else:
                paths = [expected_paths]
            if not all(exists(path) for path in paths):
                return None
            return expected_paths

    def _try_get_smdb_analysis(
        self,
        smdb: SMDB,
        target: TargetT,
        expected_path: Path,
        validate: bool,
    ) -> Path | None:
        """
        Check if SMDB already has analysis, and invalidate it if the
        output for a stage doesn't exist
        """
        if not self.analysis_type:
            return None

        analysis = target.analysis_by_type.get(self.analysis_type)

        if validate:
            ds_name = None
            if isinstance(target, Sample):
                sample = cast(Sample, target)
                sample_ids = [sample.id]
                ds_name = sample.dataset.name
            elif isinstance(target, Dataset):
                ds = cast(Dataset, target)
                sample_ids = [s.id for s in ds.get_samples()]
                ds_name = ds.name
            else:
                cohort = cast(Cohort, target)
                sample_ids = cohort.get_all_sample_ids()

            found_path = smdb.process_existing_analysis(
                sample_ids=sample_ids,
                completed_analysis=analysis,
                analysis_type=self.analysis_type.value,
                expected_output_fpath=expected_path,
                dataset_name=ds_name,
            )
        elif analysis:
            found_path = analysis.output
        else:
            found_path = None
        return found_path

    def _queue_reuse_job(
        self,
        target: TargetT,
        found_paths: Path | dict[str, Path]
    ) -> StageOutput:
        """
        Queues a [reuse] Job
        """
        attributes = {}
        if isinstance(target, Sample):
            attributes['sample'] = target.id
            attributes['dataset'] = target.dataset.name
        if isinstance(target, Dataset):
            attributes['sample'] = target.name
        return self.make_outputs(
            target=target,
            data=found_paths,
            jobs=[self.b.new_job(f'{self.name} [reuse]', attributes)]
        )
