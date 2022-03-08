"""
Stage classes
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable, cast, Union, TypeVar, Generic

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes.pipeline.analysis import AnalysisType
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes import buckets
from cpg_pipes.pipeline.cohort import Cohort
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.pair import Pair

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


StageOutputData = Union[str, Path, hb.Resource, dict[str, str], dict[str, Path], dict[str, hb.Resource]]


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
        if isinstance(data, str):
            data = Path(data)
        if isinstance(data, dict):
            data = {k: (Path(v) if isinstance(v, str) else v) for k, v in data.items()}
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
        Cast the result to str. Though exception if failed to cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, Path):
            raise ValueError(f'{res} is not a path.')
        return cast(Path, res)

    def as_resource(self, id=None) -> hb.Resource:
        """
        Cast the result to Hail Batch Resource, error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, hb.Resource):
            raise ValueError(f'{res} is not a Hail Batch Resource.')
        return cast(hb.Resource, res)

    def as_dict(self) -> dict[str, Path|hb.Resource]:
        """
        Cast the result to a dictionary, error if can't cast.
        """
        if not isinstance(self.data, dict):
            raise ValueError(f'{self.data} is not a dictionary.')
        return self.data

    def as_resource_dict(self) -> dict[str, hb.Resource]:
        """
        Cast the result to a dictionary of Hail Batch Resources, error if can't cast.
        """
        return {k: self.as_resource(id=k) for k in self.as_dict()}

    def as_path_dict(self) -> dict[str, hb.Resource]:
        """
        Cast the result to a dictionary of strings, error if can't cast.
        """
        return {k: self.as_path(id=k) for k in self.as_dict()}


StageDecorator = Callable[..., 'Stage']


class StageInput:
    """
    Represents an input for a stage run. Wraps outputs of all required upstream stages
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
            target_id = output.target.unique_id
    
            if stage_name not in self._results_by_target_by_stage:
                self._results_by_target_by_stage[stage_name] = dict()
            self._results_by_target_by_stage[stage_name][target_id] = output

    def _each(
        self, 
        fun: Callable,
        stage: StageDecorator,
    ):
        if stage.__name__ not in self._results_by_target_by_stage:
            raise ValueError(
                f'No inputs from {stage.__name__} for {self.stage.name} found '
                'after skipping targets with missing inputs. ' +
                ('Check the logs if all samples were missing inputs from previous '
                 'stages, and consider changing --first-stage'
                    if self.stage.pipe.skip_missing_input else '')
            )

        return {
            trg: fun(result)
            for trg, result 
            in self._results_by_target_by_stage.get(stage.__name__, {}).items()
        }

    def as_path_by_target(
        self, 
        stage: StageDecorator,
        id: str|None = None,
    ) -> dict[str, Path]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_path(id=id)), stage=stage)

    def as_resource_by_target(
        self, 
        stage: StageDecorator,
        id: str|None = None,
    ) -> dict[str, hb.Resource]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_resource(id=id)), stage=stage)

    def as_dict_by_target(self, stage: StageDecorator) -> dict[str, dict[str, Path]]:
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
        id: str|None = None,
    ) -> Path:
        """
        Represent as a path to a file, otherwise fail.
        `stage` can be callable, or a subclass of Stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_path(id)
    
    def as_resource(
        self, 
        target: 'Target',
        stage: StageDecorator,
        id: str|None = None,
    ) -> str:
        """
        Get Hail Batch Resource for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_resource(id)

    def as_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dictoinary of files or Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_dict()
    
    def as_path_dict(self, target: 'Target', stage: StageDecorator) -> dict[str, Path]:
        """
        Get a dictoinary of files for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_path_dict()

    def as_resource_dict(
        self, target: 'Target', stage: StageDecorator
    ) -> dict[str, Path]:
        """
        Get a dictoinary of  Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_resource_dict()

    def get_jobs(self) -> list[Job]:
        """
        Build a list of hail batch dependencies from all stages and targets
        """
        jobs = []
        for _, results_by_target in self._results_by_target_by_stage.items():
            for _, results in results_by_target.items():
                jobs.extend(results.jobs)
        return jobs


# Type variable to make sure a Stage subclass always matches the
# correspondinng Target subclass
TargetT = TypeVar('TargetT', bound=Target)

ExpectedResultT = Union[str, Path, dict[str, str], dict[str, Path], None]


class Stage(Generic[TargetT], ABC):
    """
    Abstract class for a pipeline stage.
    """
    def __init__(
        self,
        pipeline,
        name: str,
        required_stages: list[StageDecorator] | StageDecorator | None = None,
        sm_analysis_type: AnalysisType|None = None,
        skipped: bool = False,
        required: bool = True,
        assume_results_exist: bool = False,
        forced: bool = False,
    ):
        """
        :param required_stages: list of stage classes that this stage requires
        :param sm_analysis_type: if defined, will query the SMDB Analysis entries 
            of this type
        :param skipped: means that the stage is skipped and self.queue_jobs()
            won't run. The other stages if depend on it can aassume that that 
            self.expected_result() returns existing files and target.ouptut_by_stage 
            will be populated.
        :param required: means that the self.expected_output() results are 
            required for another active stage, even if the stage was skipped.
        :param assume_results_exist: for skipped but required stages, 
            the self.expected_result() output will still be checked for existence. 
            This option makes the downstream stages assume that the output exist.
        :param forced: run self.queue_jobs() even if can reuse the 
            self.expected_output().
        """
        self._name = name
        self.pipe = pipeline
        self.required_stages_classes: list[StageDecorator] = []
        if required_stages:
            if isinstance(required_stages, list):
                self.required_stages_classes.extend(required_stages)
            else:
                self.required_stages_classes.append(required_stages)

        # Populated in pipeline.run(), after we know all stages
        self.required_stages: list[Stage] = []
        
        # If analysis type is defined, it will be used to find and reuse existing
        # outputs from the SMDB
        self.analysis_type = sm_analysis_type

        # Populated with the return value of `add_to_the_pipeline()`
        self.output_by_target: dict[str, StageOutput] = dict()

        self.skipped = skipped  
        self.required = required
        self.forced = forced
        self.assume_results_exist = assume_results_exist
    
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
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `queue_jobs()`.
        
        Can be a str or a Path object, or a dictionary of str/Path objects.
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
        data: StageOutputData | None = None,
        jobs: list[Job] | None = None
    ) -> StageOutput:
        """
        Builds a StageDeps object to return from a stage's queue_jobs()
        """
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
        Performs checks like possibility to reuse existing jobs results,
        or if required dependencies are missing.
        """
        if not self.skipped:
            reusable_paths = self._try_get_reusable_paths(target)
            if reusable_paths:
                if target.forced:
                    logger.info(
                        f'{self.name}: can reuse, but forcing rerunning the stage '
                        f'for {target.unique_id}'
                    )
                    return self.queue_jobs(target, self._make_inputs())
                else:
                    logger.info(f'{self.name}: reusing results for {target.unique_id}')
                    return self._queue_reuse_job(target, reusable_paths)
            else:
                logger.info(f'{self.name}: adding jobs for {target.unique_id}')
                return self.queue_jobs(target, self._make_inputs())

        elif self.required:
            reusable_paths = self._try_get_reusable_paths(target)
            if not reusable_paths:
                if self.pipe.skip_missing_input:
                    logger.info(
                        f'Stage {self.name} is skipped and required, however '
                        f'--skip-samples-without-first-stage-input is set, so skipping '
                        f'this target ({target.unique_id}).'
                    )
                    target.active = False
                    return self.make_outputs(target=target) 
                else:
                    raise ValueError(
                        f'Stage {self.name} is required, but skipped and '
                        f'cannot reuse outputs for {target.unique_id}'
                    )
            else:
                return self.make_outputs(target=target, data=reusable_paths) 

        else:
            # Stage is not needed, returning empty outputs
            return self.make_outputs(target=target)

    def _try_get_reusable_paths(
        self, 
        target: TargetT
    ) -> dict[str, Path] | Path | None:
        """
        Returns outputs that can be reused for the stage for the target,
        or None of none can be reused
        """
        if self.pipe.dry_run:
            return None

        expected_output = self.expected_result(target)
        
        # Converting all str into Path
        expected_paths: dict[str, Path] | Path | None
        if isinstance(expected_output, dict):
            expected_paths = {k: Path(v) for k, v in expected_output.items()}
        elif isinstance(expected_output, str):
            expected_paths = Path(expected_output)
        else:
            expected_paths = expected_output

        if self.analysis_type is not None:
            # Checking the SMDB
            if not expected_output:
                raise ValueError(
                    f'expected_result() returned None, but must return str '
                    f'for a stage with analysis_type: {self.name} '
                    f'on {target.unique_id}, analysis_type={self.analysis_type}'
                )
                
            if isinstance(expected_paths, dict):
                raise ValueError(
                    f'expected_result() returns a dict, won\'t check the SMDB for '
                    f'{self.name} on {target.unique_id}'
                )
            validate = self.pipe.validate_smdb_analyses
            if self.skipped and (
                self.assume_results_exist or self.pipe.skip_missing_input
            ):
                validate = False

            # found_path would be analysis.output from DB if it exists; if it doesn't,
            # and validate_smdb_analyses=True, expected_paths file would be checked:
            found_path = self._try_get_smdb_analysis(
                target, cast(Path, expected_paths), validate=validate
            )
            return found_path
        
        elif not expected_paths:
            return None

        elif (
            not self.pipe.check_expected_outputs or
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
            if not all(buckets.file_exists(path) for path in paths):
                return None
            return expected_paths

    def _try_get_smdb_analysis(
        self, 
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

            found_path = self.pipe.db_process_existing_analysis(
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
            jobs=[self.pipe.b.new_job(f'{self.name} [reuse]', attributes)]
        )


class SampleStage(Stage[Sample], ABC):
    """
    Sample-level stage
    """
    @abstractmethod
    def expected_result(self, sample: 'Sample') -> ExpectedResultT:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """
        
    @abstractmethod
    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline) -> dict[str, StageOutput]:
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for ds_i, ds in enumerate(datasets):
            logger.info(f'{self.name}: #{ds_i}/{ds.name} {ds}')
            if not ds.get_samples():
                raise ValueError(f'No active samples are found to run in the dataset {ds.name}')
            for sample_i, sample in enumerate(ds.get_samples()):
                logger.info(f'{self.name}: #{sample_i}/{sample}')
                sample_result = self._queue_jobs_with_checks(sample)
                output_by_target[sample.unique_id] = sample_result
                logger.info('------')
            logger.info('-#-#-#-')
        return output_by_target


class PairStage(Stage, ABC):
    """
    Stage on a pair os samples
    """
    @abstractmethod
    def expected_result(self, pair: Pair) -> ExpectedResultT:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, pair: 'Pair', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline) -> dict[str, StageOutput]:
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if datasets:
            raise ValueError('No active datasets are found to run')
        for ds in datasets:
            if not ds.get_samples():
                raise ValueError(
                    f'No active samples are found to run in the dataset {ds.name}'
                )
            for s1 in ds.get_samples():
                for s2 in ds.get_samples():
                    if s1 == s2:
                        continue
                    pair = Pair(s1, s2)
                    output_by_target[pair.unique_id] = \
                        self._queue_jobs_with_checks(pair)
        return output_by_target


class DatasetStage(Stage, ABC):
    """
    Dataset-level stage
    """
    @abstractmethod
    def expected_result(self, dataset: 'Dataset') -> ExpectedResultT:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, dataset: 'Dataset', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline) -> dict[str, StageOutput]:
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for ds_i, ds in enumerate(datasets):
            logger.info(f'{self.name}: #{ds_i}/{ds.name} {ds}')
            output_by_target[ds.unique_id] = \
                self._queue_jobs_with_checks(ds)
            logger.info('-#-#-#-')
        return output_by_target


class CohortStage(Stage, ABC):    
    """
    Entire cohort level stage
    """
    @abstractmethod
    def expected_result(self, cohort: 'Cohort') -> ExpectedResultT:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, cohort: 'Cohort', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline) -> dict[str, StageOutput]:
        return {
            pipeline.cohort.unique_id: 
            self._queue_jobs_with_checks(pipeline.cohort)
        }
