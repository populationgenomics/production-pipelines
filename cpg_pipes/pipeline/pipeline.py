"""
Provides `Pipeline` class and a `@stage` decorator that allows to define pipeline
stages and plug them together.

Workflow of a stage is to add jobs to Hail Batch. Each stage acts on a `Target`, 
which can be a `Sample`, a `Dataset`, or a `Cohort` (= all input datasets combined). 
Pipeline would resolve dependencies between stages of different levels accordingly.

For examples of pipelines, see the `pipelines/` folder in the repository root.
"""
import functools
import logging
import shutil
import tempfile
import time
from abc import ABC, abstractmethod
from typing import Optional, cast, Union, Any, Callable, Type

from .analysis import AnalysisType
from .dataset import Dataset, Cohort, Sample, Pair
from .exceptions import PipelineError
from .metadata_provider import CsvMetadataProvider, MetadataProvider
from .stage import Stage, ExpectedResultT, StageInput, StageOutput
from ..cpg.smdb import SMDB, SmdbMetadataProvider
from ..cpg.storage import CPGStorageProvider
from ..hb.batch import setup_batch, Batch
from ..hb.prev_job import PrevJob
from ..storage import Namespace, StorageProvider, Path, to_path

logger = logging.getLogger(__file__)


StageDecorator = Callable[..., 'Stage']


# We record each initialised Stage subclass, so we know the default stage
# list for the case when the user doesn't pass them explicitly with set_stages()
_ALL_DEFINED_STAGES = []


def stage(
    _cls: Optional[Type[Stage]] = None, 
    *,
    sm_analysis_type: AnalysisType | None = None, 
    required_stages: list[StageDecorator] | StageDecorator | None = None,
    skipped: bool = False,
    required: bool | None = None,
    assume_results_exist: bool = False,
    forced: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with an optional argument.
    The goal is to allow cleaner defining of custom pipeline stages, without
    requiring to implement constructor. E.g.

    @stage(sm_analysis_type=AnalysisType.GVCF, required_stages=CramStage)
    class GvcfStage(SampleStage):
        def expected_result(self, sample: Sample):
            ...
        def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
            ...
    """
    def decorator_stage(cls) -> StageDecorator:
        """Implements decorator."""
        @functools.wraps(cls)
        def wrapper_stage(pipeline: 'Pipeline') -> Stage:
            """Decorator helper function."""
            return cls(
                name=cls.__name__,
                batch=pipeline.b,
                dry_run=pipeline.dry_run,
                required_stages=required_stages,
                sm_analysis_type=sm_analysis_type,
                skipped=skipped,
                required=required,
                assume_results_exist=assume_results_exist, 
                forced=forced,
                validate_smdb_analyses=pipeline.validate_smdb_analyses,
                skip_samples_with_missing_input=pipeline.skip_samples_with_missing_input,
                check_expected_outputs=pipeline.check_expected_outputs,
                check_intermediates=pipeline.check_intermediates,
                smdb=pipeline.get_db(),
                pipeline_config=pipeline.config,
            )
        # We record each initialised Stage subclass, so we know the default stage
        # list for the case when the user doesn't pass them explicitly with set_stages()
        _ALL_DEFINED_STAGES.append(wrapper_stage)
        return wrapper_stage

    if _cls is None:
        return decorator_stage
    else:
        return decorator_stage(_cls)


def skip(
    _fun: Optional[StageDecorator] = None, 
    *,
    assume_results_exist: bool = False,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Decorator on top of `@stage` that sets the `self.skipped` field to True

    @skip
    @stage
    class MyStage1(SampleStage):
        ...

    @skip
    @stage(assume_results_exist=True)
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
            s.assume_results_exist = assume_results_exist
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
        namespace: Namespace | str,
        version: str | None = None,
        storage_provider: StorageProvider | None = None,
        metadata_source: str = 'smdb',
        metadata_csv_path: str | None = None,
        stages_in_order: list[StageDecorator] | None = None,
        keep_scratch: bool = True,
        dry_run: bool = False,
        previous_batch_tsv_path: str | None = None,
        previous_batch_id: str | None = None,
        update_smdb_analyses: bool = False,
        check_smdb_seq: bool = False,
        skip_samples_with_missing_input: bool = False,
        validate_smdb_analyses: bool = False,
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
        if isinstance(namespace, str):
            namespace = Namespace(namespace)
        self.namespace = namespace
        storage_provider = storage_provider or CPGStorageProvider()
        self.storage_provider = storage_provider

        self.cohort = Cohort(
            self.name,
            analysis_dataset_name=analysis_dataset,
            namespace=namespace,
            storage_provider=storage_provider
        )

        self.check_intermediates = check_intermediates and not dry_run
        self.check_expected_outputs = check_expected_outputs and not dry_run
        self.skip_samples_with_missing_input = skip_samples_with_missing_input
        self.first_stage = first_stage
        self.last_stage = last_stage
        self.validate_smdb_analyses = validate_smdb_analyses

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

        self.prev_batch_jobs = dict()
        if previous_batch_tsv_path is not None:
            assert previous_batch_id is not None
            self.prev_batch_jobs = PrevJob.parse(
                previous_batch_tsv_path,
                previous_batch_id,
                tmp_bucket=self.tmp_bucket,
                keep_scratch=keep_scratch,
            )

        self.config = config or {}

        if version:
            description += f' {version}'
        if input_datasets:
            description += ': ' + ', '.join(input_datasets)
        self.b: Batch = setup_batch(
            description=description, 
            tmp_bucket=self.tmp_bucket,
            keep_scratch=self.keep_scratch,
            billing_project=self.cohort.analysis_dataset.stack,
        )

        self._db = None
        if input_datasets:
            metadata_provider: MetadataProvider
            if metadata_source == 'smdb':
                self._db = SMDB(
                    self.cohort.analysis_dataset.name,
                    do_update_analyses=update_smdb_analyses,
                )
                metadata_provider = SmdbMetadataProvider(self._db)
            else:
                assert metadata_source == 'csv'
                if not metadata_csv_path:
                    raise PipelineError(
                        '--metadata-source is "csv", --metadata-csv path must be '
                        'provided'
                    )
                with to_path(metadata_csv_path).open() as fp:
                    metadata_provider = CsvMetadataProvider(fp=fp)

            metadata_provider.populate_cohort(
                cohort=self.cohort,
                dataset_names=input_datasets,
                skip_samples=skip_samples,
                only_samples=only_samples,
                ped_files=ped_files,
                do_check_seq_existence=check_smdb_seq,
            )

        if force_samples:
            for s in self.cohort.get_all_samples():
                if s.id in force_samples:
                    logger.info(
                        f'Force rerunning sample {s.id} even if its outputs exist'
                    )
                    s.forced = True

        self._stages_dict: dict[str, Stage] = dict()
        self._stages_in_order: list[StageDecorator] = stages_in_order or _ALL_DEFINED_STAGES

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
        # TODO: use TopologicalSorter
        # import graphlib
        # for stage in graphlib.TopologicalSorter(stages).static_order():
        # https://github.com/populationgenomics/analysis-runner/pull/328/files
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
                # breakpoint()
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

    @property
    def db(self) -> SMDB:
        """
        Get read-only sample-metadata DB object.
        """
        if self._db is None:
            raise PipelineError('SMDB is not initialised')
        return cast('SMDB', self._db)

    def get_db(self) -> SMDB | None:
        """
        Like .db property, but returns None if db is not initialised.
        """
        return self._db

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
        return self.cohort.get_all_samples(only_active=only_active)

    def get_all_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Thin wrapper around corresponding Cohort method.
        """
        return self.cohort.get_all_sample_ids(only_active=only_active)


class SampleStage(Stage[Sample], ABC):
    """
    Sample-level stage.
    """

    @abstractmethod
    def expected_result(self, sample: Sample) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Pplug in stage into the pipeline.
        """
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for ds_i, ds in enumerate(datasets):
            logger.info(f'{self.name}: #{ds_i} {ds}')
            if not ds.get_samples():
                raise ValueError(
                    f'No active samples are found to run in the dataset {ds.name}')
            for sample_i, sample in enumerate(ds.get_samples()):
                logger.info(f'{self.name}: #{sample_i}/{sample}')
                sample_result = self._queue_jobs_with_checks(sample)
                output_by_target[sample.target_id] = sample_result
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
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, pair: Pair, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Pplug in stage into the pipeline.
        """
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
                    output_by_target[pair.target_id] = \
                        self._queue_jobs_with_checks(pair)
        return output_by_target


class DatasetStage(Stage, ABC):
    """
    Dataset-level stage
    """

    @abstractmethod
    def expected_result(self, dataset: Dataset) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Pplug in stage into the pipeline.
        """
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for ds_i, ds in enumerate(datasets):
            logger.info(f'{self.name}: #{ds_i}/{ds.name} {ds}')
            output_by_target[ds.target_id] = \
                self._queue_jobs_with_checks(ds)
            logger.info('-#-#-#-')
        return output_by_target


class CohortStage(Stage, ABC):
    """
    Entire cohort level stage
    """

    @abstractmethod
    def expected_result(self, cohort: Cohort) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Override to plug in stage into the pipeline.
        """
        return {
            pipeline.cohort.target_id:
                self._queue_jobs_with_checks(pipeline.cohort)
        }
