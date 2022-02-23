"""
Provides "Pipeline" class that allows plugging multiple "stages" together
by resolving dependencies through the sample-metadata database, or by checking
objects on buckets directly. 

Each stage adds jobs to Hail Batch. Each stage acts on "target", which can be a 
sample, a project, or an entire cohort. Pipeline would resolve dependencies between
stages of different levels accordingly.

Basic example:

@stage(analysis_type=AnalysisType.CRAM)
class CramStage(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        expected_path = self.expected_result(pipe)
        job = align.bwa(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, data=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.GVCF, requires_stages=CramStage)
class GvcfStage(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = inputs.as_path(target=sample, stage=CramStage)
        expected_path = self.expected_result(pipe)
        job = haplotype_caller.produce_gvcf(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, data=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.JOINT_CALLING)
class JointCallingStage(CohortStage):
    def queue_jobs(self, pipeline: Pipeline, inputs: StageInput) -> StageOutput:
        gvcf_by_sid = inputs.as_path_by_target(stage=GvcfStage)
        expected_path = self.expected_result(pipe)
        job = make_joint_genotyping_jobs(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(pipe, data=expected_path, jobs=[job])

@click.command()
@pipeline_click_options
def main(**kwargs):
    p = Pipeline(
        name='my_joint_calling_pipeline',
        title='My joint calling pipeline',
        stages_in_order=[CramStage, GvcfStage, JointCallingStage],
        **kwargs
    )
    p.submit_batches()

For more usage examples, see the "pipelines" folder in the root of this repository.
"""
import functools
import logging
import shutil
import tempfile
from typing import List, Dict, Optional, Tuple, cast, Union, Any, Callable, Type

from cpg_pipes import buckets
from cpg_pipes.hb.batch import setup_batch, Batch
from cpg_pipes.hb.prev_job import PrevJob
from cpg_pipes.namespace import Namespace
from cpg_pipes.pipeline.cohort import Cohort
from cpg_pipes.pipeline.project import Project
from cpg_pipes.pipeline.stage import Stage
from cpg_pipes.smdb.types import AnalysisType
from cpg_pipes.smdb.smdb import SMDB


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


StageDecorator = Callable[..., 'Stage']


# We record each initialised Stage subclass, so we know the default stage
# list for the case when the user doesn't pass them explicitly with set_stages()
_defined_stages = []


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
        # We record each initialised Stage subclass, so we know the default stage
        # list for the case when the user doesn't pass them explicitly with set_stages()
        global _defined_stages
        _defined_stages.append(wrapper_stage)
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


def run_pipeline(dry_run: bool = False, **kwargs) -> 'Pipeline':
    """
    Create and submit a pipeline to Hail Batch.
    """
    pipeline = Pipeline(**kwargs)
    pipeline.submit_batch(dry_run=dry_run)
    return pipeline


class PipelineException(Exception):
    pass


class Pipeline:
    """
    Represents a Pipeline, and incapulates a Hail Batch object,
    stages, and a cohort of projects with samples.
    """
    def __init__(
        self,
        analysis_project: str,
        name: str,
        title: str,
        output_version: str,
        namespace: Union[Namespace, str],
        stages_in_order: Optional[List[StageDecorator]] = None,
        keep_scratch: bool = True,
        dry_run: bool = False,
        previous_batch_tsv_path: Optional[str] = None,
        previous_batch_id: Optional[str] = None,
        update_smdb_analyses: bool = False,
        check_smdb_seq_existence: bool = False,
        skip_samples_without_first_stage_input: bool = False,
        validate_smdb_analyses: bool = False,
        check_intermediate_existence: bool = True,
        check_job_expected_outputs_existence: bool = True,
        first_stage: Optional[str] = None,
        last_stage: Optional[str] = None,
        config: Optional[Dict] = None,
        input_projects: Optional[List[str]] = None,
        source_tag: Optional[str] = None,
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        force_samples: Optional[List[str]] = None,
        ped_files: Optional[List[str]] = None,
        local_dir: Optional[str] = None,
    ):
        if stages_in_order and not input_projects:
            raise ValueError(
                'Projects must be populated before adding stages. '
                'Provide `input_projects`, or omit `stages_in_order` and call '
                'pipeline.set_stages(stages_in_order) later.'
            )

        super().__init__()
        if isinstance(namespace, str):
            namespace = Namespace(namespace)
        self.analysis_project = Project(
            name=analysis_project,
            namespace=namespace,
            pipeline=self,
        )
        self.name = name
        self.output_version = output_version
        self.namespace = namespace
        self.check_intermediate_existence = check_intermediate_existence and not dry_run
        self.check_job_expected_outputs_existence = check_job_expected_outputs_existence and not dry_run
        self.skip_samples_without_first_stage_input = skip_samples_without_first_stage_input
        self.first_stage = first_stage
        self.last_stage = last_stage
        self.validate_smdb_analyses = validate_smdb_analyses

        if namespace == Namespace.TMP:
            tmp_suf = 'test-tmp'
            analysis_suf = 'test-tmp/analysis'
            web_suf = 'test-tmp/web'
            self.output_suf = 'test-tmp'
            self.proj_output_suf = 'test'
        elif namespace == Namespace.TEST:
            tmp_suf = 'test-tmp'
            analysis_suf = 'test-analysis'
            web_suf = 'test-web'
            self.output_suf = 'test'
            self.proj_output_suf = 'test'
        else:
            tmp_suf = 'main-tmp'
            analysis_suf = 'main-analysis'
            web_suf = 'main-web'
            self.output_suf = 'main'
            self.proj_output_suf = 'main'
        
        path_ptrn = (
            f'gs://cpg-{self.analysis_project.stack}-{{suffix}}/'
            f'{self.name}/'
            f'{self.output_version}'
        )
        self.tmp_bucket = path_ptrn.format(suffix=tmp_suf)
        self.analysis_bucket = path_ptrn.format(suffix=analysis_suf)
        self.web_bucket = path_ptrn.format(suffix=web_suf)
        self.web_url = (
            f'https://{self.namespace.value}-web.populationgenomics.org.au/'
            f'{self.analysis_project.stack}/'
            f'{self.name}/'
            f'{self.output_version}'
        )
        self.keep_scratch = keep_scratch
        self.dry_run: bool = dry_run

        if local_dir:
            self.local_dir = local_dir
            self.local_tmp_dir = tempfile.mkdtemp(dir=local_dir)
        else:
            self.local_dir = self.local_tmp_dir = tempfile.mkdtemp()

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

        self.b: Batch = setup_batch(
            title=title, 
            tmp_bucket=self.tmp_bucket,
            keep_scratch=self.keep_scratch,
            billing_project=self.analysis_project.stack,
        )

        self.cohort = Cohort(self.name, pipeline=self)
        self._db = None
        if input_projects:
            self._db = SMDB(
                self.analysis_project.name,
                do_update_analyses=update_smdb_analyses,
                do_check_seq_existence=check_smdb_seq_existence,
            )
            self.cohort.populate(
                smdb=self._db,
                input_projects=input_projects,
                local_tmp_dir=self.local_dir,
                source_tag=source_tag,
                skip_samples=skip_samples,
                only_samples=only_samples,
                ped_files=ped_files,
                forced_samples=force_samples,
            )

        self._stages_dict: Dict[str, Stage] = dict()
        if stages_in_order:
            self.set_stages(stages_in_order)
        else:
            self.set_stages(_defined_stages)

    def submit_batch(
        self, 
        dry_run: Optional[bool] = None,
        wait: Optional[bool] = False,
    ) -> Optional[Any]:
        """
        Submits Hail Batch jobs.
        """
        if dry_run is None:
            dry_run = self.dry_run
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

    def _validate_first_last_stage(self) -> Tuple[Optional[int], Optional[int]]:
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

    def set_stages(self, stages_classes: List[StageDecorator]):
        """
        Iterate over stages and call add_to_the_pipeline() on each.
        Effectively creates all Hail Batch jobs through Stage.queue_jobs().
        """
        # Initializing stage objects
        stages = [cls(self) for cls in stages_classes]
        for stage in stages:
            if stage.name in self._stages_dict:
                raise ValueError(
                    f'Stage {stage.name} is already defined. Check your '
                    f'list for duplicates: {", ".join(s.name for s in stages)}'
                )
            self._stages_dict[stage.name] = stage

        first_stage_num, last_stage_num = self._validate_first_last_stage()

        # First round - checking which stages we require, even if they are skipped.
        # If there are required stages that are not defined, we defined them
        # into `additional_stages` as `skipped`.
        # TODO: use TopologicalSorter
        # import graphlib
        # for stage in graphlib.TopologicalSorter(stages).static_order():
        # https://github.com/populationgenomics/analysis-runner/pull/328/files
        additional_stages_dict: Dict[str, Stage] = dict()
        for i, (stage_name, stage) in enumerate(self._stages_dict.items()):
            if first_stage_num is not None and i < first_stage_num:
                stage.skipped = True
                stage.required = False
                logger.info(f'Skipping stage {stage_name}')
                continue

            if last_stage_num is not None and i > last_stage_num:
                stage.skipped = True
                stage.required = False
                continue

            for reqcls in stage.required_stages_classes:
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
                        f'but the output will be required for the stage {stage.name}'
                    )
                stage.required_stages.append(reqstage)
        
        if additional_stages_dict:
            logger.info(
                f'Additional required stages added as "skipped": '
                f'{additional_stages_dict.keys()}'
            )
            for name, stage in self._stages_dict.items():
                additional_stages_dict[name] = stage    
            self._stages_dict = additional_stages_dict

        logger.info(f'Setting stages: {", ".join(s.name for s in stages if not s.skipped)}')

        # Second round - actually adding jobs from the stages.
        for i, stage in enumerate(self._stages_dict.values()):
            if not stage.skipped:
                logger.info(f'*' * 60)
                logger.info(f'Stage {stage.name}')

            if stage.required:
                logger.info(f'Adding jobs for stage {stage.name}')
                stage.output_by_target = stage.add_to_the_pipeline(self)

            if not stage.skipped:
                logger.info(f'')
                if last_stage_num and i >= last_stage_num:
                    logger.info(f'Last stage is {stage.name}, stopping here')
                    break

    def can_reuse(self, fpath: Optional[str]) -> bool:
        """
        Checks if the fpath exists, 
        but always returns False if not check_intermediate_existence
        """
        return buckets.can_reuse(fpath, overwrite=not self.check_intermediate_existence)

    def db_process_existing_analysis(self, *args, **kwargs):
        if self._db is None:
            raise PipelineException('SMDB is not initialised')
        
        return self._db.process_existing_analysis(*args, **kwargs)

    @property
    def db(self) -> SMDB:
        if self._db is None:
            raise PipelineException('SMDB is not initialised')
        return cast('SMDB', self._db)

    def get_db(self) -> Optional[SMDB]:
        """
        Like .db property, but returns None if db is not initiazlied
        """
        return self._db
