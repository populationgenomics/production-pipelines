"""
Implements a Pipeline that is capable of plugging multiple "stages" together,
by resolving dependencies through the sample-metadata database or by checking
buckets directly. Each stage adds jobs to Hail Batch. Each stage acts on "target",
which can be a sample, a project, or an entire cohort. Pipeline would resolve 
dependencies between stages of different levels accordingly.

For usage examples, see the "pipelines" folder in the root of this repository.

Basic example is also provided here:

@stage(analysis_type=AnalysisType.CRAM)
class CramStage(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInputs) -> StageOutputs:
        expected_path = self.get_expected_output(pipe)
        job = align.bwa(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, data=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.GVCF, requires_stages=CramStage)
class GvcfStage(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInputs) -> StageOutputs:
        cram_path = inputs.as_path(target=sample, stage=CramStage)
        expected_path = self.get_expected_output(pipe)
        job = haplotype_caller.produce_gvcf(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, data=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.JOINT_CALLING)
class JointCallingStage(CohortStage):
    def queue_jobs(self, pipeline: Pipeline, inputs: StageInputs) -> StageOutputs:
        gvcf_by_sid = inputs.as_path_by_target(stage=GvcfStage)
        expected_path = self.get_expected_output(pipe)
        job = make_joint_genotyping_jobs(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(pipe, data=expected_path, jobs=[job])

@click.command()
@pipeline_click_options
def main(**kwargs):
    run_pipeline(
        name='my_joint_calling_pipeline',
        title='My joint calling pipeline',
        stages_in_order=[CramStage, GvcfStage, JointCallingStage],
        **kwargs
    )
"""

import functools
import shutil
import tempfile
from dataclasses import dataclass
import logging
from enum import Enum
from os.path import join
from typing import List, Dict, Optional, Tuple, Callable, cast, Type, Union, TypeVar, \
    Generic
from abc import ABC, abstractmethod

import click
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import utils
from cpg_pipes.utils import Namespace, AnalysisType
from cpg_pipes.smdb import SMDB, Analysis, Sequence
from cpg_pipes.hailbatch import AlignmentInput, PrevJob, get_hail_bucket, setup_batch

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def pipeline_click_options(function: Callable) -> Callable:
    """
    Decorator to use with click when writing a script that implements a pipeline.
    For example:

    @click.command()
    @pipeline_click_options
    @click.argument('--custom-argument')
    def main():
        pass
    """
    options = [
        click.option(
            '-n',
            '--namespace',
            'namespace',
            type=click.Choice([n.lower() for n in Namespace.__members__]),
            callback=lambda c, p, v: getattr(Namespace, v.upper()) if v else None,
            help='The bucket namespace to write the results to',
        ),
        click.option(
            '--analysis-project',
            'analysis_project',
            default='seqr',
            help='SM project name to write the intermediate/joint-calling analysis '
                 'entries to',
        ),
        click.option(
            '--input-project',
            'input_projects',
            multiple=True,
            required=True,
            help='Only read samples that belong to the project(s). Can be set multiple '
                 'times.',
        ),
        click.option(
            '--ped-file',
            'ped_files',
            multiple=True,
        ),
        click.option(
            '--first-stage',
            'first_stage',
            help='Only pick results from the previous stages if they exist. '
            'If not, skip such samples',
        ),
        click.option(
            '--last-stage',
            'last_stage',
            help='Finish the pipeline after this stage',
        ),
        click.option(
            '--skip-sample',
            '-S',
            'skip_samples',
            multiple=True,
            help='Don\'t process specified samples. Can be set multiple times.',
        ),
        click.option(
            '--only-sample',
            '-s',
            'only_samples',
            multiple=True,
            help='Only take these samples (can be set multiple times)',
        ),
        click.option(
            '--force-sample',
            'force_samples',
            multiple=True,
            help='Force reprocessing these samples. Can be set multiple times.',
        ),
        click.option(
            '--output-version',
            'output_version',
            type=str,
            default='v0',
            help='Suffix the outputs with this version tag. Useful for testing',
        ),
        click.option(
            '--keep-scratch/--remove-scratch', 
            'keep_scratch', 
            default=True,
            is_flag=True,
        ),
        click.option('--dry-run', 'dry_run', is_flag=True),
        click.option(
            '--check-smdb-seq-existence/--no-check-smdb-seq-existence',
            'check_smdb_seq_existence',
            default=False,
            is_flag=True,
            help='Check that files in sequence.meta exist'
        ),
        click.option(
            '--skip-samples-without-seq-input',
            'skip_samples_without_seq_input',
            default=False,
            is_flag=True,
            help='If sequence.meta files for a sample don\'t exist, remove this sample '
                 'instead of failing'
        ),
        click.option(
            '--check-intermediate-existence/--no-check-intermediate-existence',
            'check_intermediate_existence',
            default=True,
            is_flag=True,
            help='Before running a job, check for an intermediate output before '
                 'submitting it, and if it exists on a bucket, submit a [reuse] job '
                 'instead. Works well with '
                 '--previous-batch-tsv/--previous-batch-id options.',
        ),
        click.option(
            '--update-smdb-analyses/--no-update-smdb-analyses',
            'update_smdb_analyses',
            is_flag=True,
            default=True,
            help='Create analysis entries for queued/running/completed jobs'
        ),
        click.option(
            '--validate-smdb-analyses/--no-validate-smdb-analyses',
            'validate_smdb_analyses',
            is_flag=True,
            default=False,
            help='Validate existing analysis entries by checking if a.output exists on '
                 'the bucket. Set the analysis entry to "failure" if output doesn\'t '
                 'exist'
        ),
        click.option(
            '--previous-batch-tsv',
            'previous_batch_tsv_path',
            help='A list of previous successful attempts from another batch, dumped '
                 'from from the Batch database (the "jobs" table joined on '
                 '"job_attributes"). If the intermediate output for a job exists in '
                 'a previous attempt, it will be passed forward, and a [reuse] job will '
                 'be submitted.'
        ),
        click.option(
            '--previous-batch-id',
            'previous_batch_id',
            help='6-letter ID of the previous successful batch (corresponds to the '
                 'directory name in the batch logs. e.g. feb0e9 in '
                 'gs://cpg-seqr-main-tmp/hail/batch/feb0e9'
        )
    ]
    for opt in options:
        function = opt(function)
    return function


StageDecorator = Callable[..., 'Stage']


def stage(
    _cls: Optional[Type['Stage']] = None, 
    *,
    sm_analysis_type: Optional[AnalysisType] = None, 
    requires_stages: Optional[Union[List[StageDecorator], StageDecorator]] = None,
) -> Union[StageDecorator, Callable[..., StageDecorator]]:
    """
    Implements a standard class decorator pattern with an optional argument.
    The goal is to allow cleaner defining of custom pipeline stages, without
    requiring to implement constructor. E.g.

    @stage(sm_analysis_type=AnalysisType.GVCF, requires_stages=CramStage)
    class GvcfStage(SampleStage):
        def get_expected_output(self, sample: Sample):
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
            )
        return wrapper_stage

    if _cls is None:
        return decorator_stage
    else:
        return decorator_stage(_cls)
    

StageOutputData = Union[str, hb.Resource, Dict[str, str], Dict[str, hb.Resource]]


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
        jobs: Optional[List[Job]] = None,
    ):
        self.data = data
        self.stage = stage
        self.target = target
        self.jobs: List[Job] = jobs or []

    def exist(self) -> bool:
        """
        Checks file existence. For Hail Batch Resources, just returns True.
        """
        if isinstance(self.data, str):
            return utils.file_exists(self.data)
        elif isinstance(self.data, dict):
            return all(utils.file_exists(v) for v in self.data.values())
        return True
    
    def __repr__(self) -> str:
        return f'result {self.data} for target {self.target}, stage {self.stage}'

    def as_path_or_resource(self, id=None) -> Union[str, hb.Resource]:
        """
        Cast the result to Union[str, hb.Resource], error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        if id is not None:
            if not isinstance(self.data, dict):
                raise ValueError(
                    f'{self} is not a dictionary, can\'t get "{id}".'
                )
            return cast(dict, self.data)[id]
        
        if isinstance(self.data, dict):
            res = cast(dict, self.data)
            if len(res.values()) > 1:
                raise ValueError(
                    f'{self} is a dictionary with more than 1 element, '
                    f'please set the `id` parameter'
                )
            return list(res.values())[0]

        return self.data        

    def as_path(self, id=None) -> str:
        """
        Cast the result to str, error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, str):
            raise ValueError(f'{self} is not a path.')
        return cast(str, res)

    def as_resource(self, id=None) -> hb.Resource:
        """
        Cast the result to Hail Batch Resource, error if can't cast.
        `id` is used to extract the value when the result is a dictionary.
        """
        res = self.as_path_or_resource(id)
        if not isinstance(res, hb.Resource):
            raise ValueError(f'{self} is not a Hail Batch Resource.')
        return cast(hb.Resource, res)

    def as_dict(self) -> Dict[str, Union[str, hb.Resource]]:
        """
        Cast the result to a dictionary, error if can't cast.
        """
        if not isinstance(self.data, dict):
            raise ValueError(f'{self} is not a dictionary.')
        return self.data
    
    def as_resource_dict(self) -> Dict[str, hb.Resource]:
        """
        Cast the result to a dictionary of Hail Batch Resources, error if can't cast.
        """
        return {k: self.as_resource(id=k) for k in self.as_dict()}

    def as_path_dict(self) -> Dict[str, hb.Resource]:
        """
        Cast the result to a dictionary of strings, error if can't cast.
        """
        return {k: self.as_path(id=k) for k in self.as_dict()}


class StageInput:
    """
    Represents an input for a stage run. Wraps outputs of all required upstream stages
    for corresponding targets (e.g. all GVCFs from a HaploytypeCallerStage
    for a JointCallingStage, along with Hail Batch jobs).

    An object of this class is passed to the public `queue_jobs` method of a Stage, 
    and can be used to query dependency files and jobs.
    """
    def __init__(self):
        self._results_by_target_by_stage: Dict[str, Dict[str, StageOutput]] = {}
        self._jobs: List[Job] = []

    def add_other_stage_output(self, output: StageOutput):
        """
        Add output from another stage run
        """
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
        return {
            trg: fun(result)
            for trg, result 
            in self._results_by_target_by_stage[stage.__name__].items()
        }

    def as_path_by_target(
        self, 
        stage: StageDecorator,
        id: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_path(id=id)), stage=stage)

    def as_resource_by_target(
        self, 
        stage: StageDecorator,
        id: Optional[str] = None,
    ) -> Dict[str, hb.Resource]:
        """
        Get a single file path result, indexed by target for a specific stage
        """
        return self._each(fun=(lambda r: r.as_resource(id=id)), stage=stage)

    def as_dict_by_target(self, stage: StageDecorator) -> Dict[str, Dict[str, str]]:
        """
        Get as a dictoinary of files/resources for a specific stage, indexed by target
        """
        return self._each(fun=(lambda r: r.as_dict(id=id)), stage=stage)

    def as_resource_dict_by_target(
        self, 
        stage: StageDecorator,
    ) -> Dict[str, Dict[str, hb.Resource]]:
        """
        Get a dictoinary of resources for a specific stage, and indexed by target
        """
        return self._each(fun=(lambda r: r.as_resource_dict(id=id)), stage=stage)

    def as_path_dict_by_target(
        self, 
        stage: StageDecorator,
    ) -> Dict[str, Dict[str, str]]:
        """
        Get a dictoinary of paths for a specific stage, and indexed by target
        """
        return self._each(fun=(lambda r: r.as_path_dict(id=id)), stage=stage)

    def as_path(
        self, 
        target: 'Target',
        stage: StageDecorator,
        id: Optional[str] = None,
    ) -> str:
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
        id: Optional[str] = None,
    ) -> str:
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_resource(id)
    
    def as_dict(self, target: 'Target', stage: StageDecorator) -> Dict[str, str]:
        """
        Get a dictoinary of files or Resources for a specific target and stage
        """
        res = self._results_by_target_by_stage[stage.__name__][target.unique_id]
        return res.as_dict()

    def get_jobs(self) -> List[Job]:
        """
        Build a list of hail batch dependencies from all stages and targets
        """
        jobs = []
        for _, results_by_target in self._results_by_target_by_stage.items():
            for _, results in results_by_target.items():
                jobs.extend(results.jobs)
        return jobs


class Sex(Enum):
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2


@dataclass
class PedigreeInfo:
    fam_id: str
    dad: Optional['Sample']
    mom: Optional['Sample']
    sex: Sex
    phenotype: str


class Target:
    """
    Stage target: e.g. sample, a dataset/project, cohort/pipeline.
    """
    def __init__(self):
        # From SMDB Analysis entries:
        self.analysis_by_type: Dict[AnalysisType, Analysis] = dict()
        # Whether to process even if outputs exist:
        self.forced = False

    @abstractmethod
    @property
    def unique_id(self) -> str:
        """
        ID should be unique across target of all levels.
        """
        pass
 

@dataclass(init=False)
class Sample(Target):
    """
    Corresponds to one Sample entry in the SMDB.
    """
    id: str
    external_id: str
    project: 'Project'
    alignment_input: Optional[AlignmentInput] = None
    good: bool = True
    seq_info: Optional[Sequence] = None
    pedigree: Optional[PedigreeInfo] = None

    def __init__(
        self, 
        id: str, 
        external_id: str, 
        project: 'Project',
    ):
        super().__init__()
        self.id = id
        self.external_id = external_id
        self.project = project

    def get_ped_dict(self, use_ext_id: bool = False) -> Dict:
        """
        Returns a dictionary of pedigree fields for this sample
        """
        def _get_id(_s: Optional[Sample]):
            if _s is None:
                return '0'
            elif use_ext_id:
                return _s.external_id
            else:
                return _s.id

        if self.pedigree:
            return {
                'Family.ID': self.pedigree.fam_id,
                'Individual.ID': _get_id(self),
                'Father.ID': _get_id(self.pedigree.dad),
                'Mother.ID': _get_id(self.pedigree.mom),
                'Sex': str(self.pedigree.sex.value),
                'Phenotype': self.pedigree.phenotype,
            }
        else:
            return {
                'Family.ID': _get_id(self),
                'Individual.ID': _get_id(self),
                'Father.ID': '0',
                'Mother.ID': '0',
                'Sex': '0',
                'Phenotype': '0',
            }

    @property
    def unique_id(self) -> str:
        return self.id


class Pair(Target):
    """
    Pair of samples
    """
    @property
    def unique_id(self) -> str:
        return f'{self.s1.unique_id}:{self.s2.unique_id}'

    def __init__(self, s1: 'Sample', s2: 'Sample'):
        super().__init__()
        self.s1 = s1
        self.s2 = s2


class Project(Target):
    """
    Represents a CPG stack (aka dataset), or a project in the sample metadata database.
    """
    def get_bucket(self):
        """
        The primary project bucket (-main or -test) 
        """
        return f'gs://cpg-{self.stack}-{self.pipeline.output_suf}'

    def get_tmp_bucket(self):
        """
        The tmp bucket (-main-tmp or -test-tmp)
        """
        return (
            f'gs://cpg-{self.stack}-{self.pipeline.output_suf}-tmp/'
            f'{self.pipeline.name}/'
            f'{self.pipeline.output_version}'
        )
    
    def __repr__(self):
        return self.name

    def __init__(
        self, 
        pipeline: 'Pipeline',
        name: str, 
        namespace: Namespace
    ):
        """
        Has "name" and "stack".
        "name" is in the SMDB terms: e.g. can be "seqr", "seqr-test".
        "stack" is in the dataset terms, so can be only e.g. "seqr", but not "seqr-test".
        """
        super().__init__()
        self.pipeline = pipeline
        self.samples: List[Sample] = []
        self.is_test = namespace != Namespace.MAIN

        if name.endswith('-test'):
            self.is_test = True
            self.stack = name[:-len('-test')]
        else:
            self.stack = name
        
        self.name = self.stack
        if self.is_test:
            self.name = self.stack + '-test'

    @property
    def unique_id(self) -> str:
        return self.name


# Type variable to make sure a Stage subclass always matches the
# correspondinng Target subclass
TargetT = TypeVar('TargetT', bound=Target)


class Stage(Generic[TargetT], ABC):
    """
    Abstract class for a pipeline stage.
    """
    def __init__(
        self,
        pipeline: 'Pipeline',
        name: str,
        requires_stages: Optional[Union[List[StageDecorator], StageDecorator]] = None,
        sm_analysis_type: Optional[AnalysisType] = None,
    ):
        self._name = name
        self.pipe = pipeline
        self.required_stages_classes: List[StageDecorator] = []
        if requires_stages:
            if isinstance(requires_stages, list):
                self.required_stages_classes.extend(requires_stages)
            else:
                self.required_stages_classes.append(requires_stages)

        # Populated in pipeline.run(), after we know all stages
        self.required_stages: List[Stage] = []
        
        # If analysis type is defined, it will be used to find and reuse existing
        # outputs from the SMDB
        self.analysis_type = sm_analysis_type

        # Populated with the return value of `add_to_the_pipeline()`
        self.output_by_target: Dict[str, StageOutput] = dict()

        # self.active=True means that the stage wasn't skipped and jobs were added
        # into the pipeline, and we shoule expect target.ouptut_by_stage 
        # to be populated. Otherwise, self.get_expected_output() should work.
        self.skipped = False  
        # self.required=True means that is required for another active stage,
        # even if it was skipped
        self.required = True
        # Rerun the stage even if can reuse the results
        self.forced = False
    
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
    def get_expected_output(self, target: TargetT) -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `queue_jobs()`.
        """
    
    @abstractmethod
    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> Dict[str, StageOutput]:
        """
        Calls `output = pipeline.add_for_target(target)` on each target, 
        which itself calls `output = queue_jobs(target, input)`, making sure to
        construct the correct `input`.

        Returns a dictionary of StageOutput, indexed by target unique_id.
        """

    def make_outputs(
        self, 
        target: TargetT,
        data: Optional[StageOutputData] = None,
        jobs: Optional[List[Job]] = None
    ) -> StageOutput:
        """
        Builds a StageDeps object to return from a stage's queue_jobs()
        """
        return StageOutput(stage=self, target=target, data=data, jobs=jobs)

    def _make_inputs(self) -> StageInput:
        """
        Collects outputs from all dependencies and create input for this stage
        """
        inputs = StageInput()
        for prev_stage in self.required_stages:
            for _, stage_output in prev_stage.output_by_target.items():
                inputs.add_other_stage_output(stage_output)
        return inputs
    
    def _queue_jobs_with_checks(self, target: TargetT) -> StageOutput:
        """
        Constructs `inputs` and calls the public `output = queue_jobs(target, input)`.
        Performs checks th like possibility to reuse existing jobs results,
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
                raise ValueError(
                    f'Stage {self.name} is required, but skipped and '
                    f'cannot reuse outputs for {target.unique_id}'
                )
            return self.make_outputs(target=target, data=reusable_paths) 

        else:
            # Stage is not needed, returning empty outputs
            return self.make_outputs(target=target)

    def _try_get_reusable_paths(
        self, 
        target: TargetT
    ) -> Optional[Union[Dict[str, str], str]]:
        """
        Returns outputs that can be reused for the stage for the target,
        or None of none can be reused
        """
        expected_output = self.get_expected_output(target)

        if self.analysis_type is not None:
            if not expected_output:
                raise ValueError(
                    f'get_expected_output() returned None, but must return str '
                    f'for a stage with analysis_type: {self.name} '
                    f'on {target.unique_id}, analysis_type={self.analysis_type}'
                )
                
            if isinstance(expected_output, dict):
                raise ValueError(
                    f'get_expected_output() returns a dict, won\'t check the SMDB for '
                    f'{self.name} on {target.unique_id}'
                )
            found_path = self._try_get_smdb_analysis(target, cast(str, expected_output))
            if found_path and not self.pipe.validate_smdb_analyses:
                return found_path

        if expected_output and self.pipe.check_intermediate_existence:
            if isinstance(expected_output, dict):
                paths = list(expected_output.values())
            else:
                paths = [expected_output]
            if all(utils.file_exists(path) for path in paths):
                return expected_output
        return None

    def _try_get_smdb_analysis(self, target: TargetT, expected_path: str) -> Optional[str]:
        """
        Check if SMDB already has analysis, and invalidate it if the
        output for a stage doesn't exist
        """
        if not self.analysis_type:
            return None
        analysis = target.analysis_by_type.get(self.analysis_type)
        if not analysis:
            return None

        if self.pipe.validate_smdb_analyses:
            if isinstance(target, Sample):
                sample = cast(Sample, target)
                sample_ids = [sample.id]
            elif isinstance(target, Project):
                project = cast(Project, target)
                sample_ids = [s.id for s in project.samples]
            else:
                pipe = cast(Pipeline, target)
                sample_ids = pipe.get_all_sample_ids()

            found_path = self.pipe.db.process_existing_analysis(
                sample_ids=sample_ids,
                completed_analysis=analysis,
                analysis_type=self.analysis_type.value,
                expected_output_fpath=expected_path,
            )
        else:
            found_path = analysis.output
        return found_path

    def _queue_reuse_job(
        self,
        target: TargetT,
        found_paths: Union[str, Dict[str, str]]
    ) -> StageOutput:
        """
        Queues a [reuse] Job
        """
        attributes = {}
        if isinstance(target, Sample):
            attributes['sample'] = target.id
            attributes['project'] = target.project.name
        if isinstance(target, Project):
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
    def get_expected_output(
        self, 
        sample: 'Sample'
    ) -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """
        
    @abstractmethod
    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> Dict[str, StageOutput]:
        output_by_target = dict()
        for project in pipeline.projects:
            for sample in project.samples:
                sample_result = self._queue_jobs_with_checks(sample)
                output_by_target[sample.unique_id] = sample_result
        return output_by_target


class PairStage(Stage, ABC):
    """
    Stage on a pair os samples
    """
    @abstractmethod
    def get_expected_output(
        self, 
        pair: Pair
    ) -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, pair: 'Pair', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> Dict[str, StageOutput]:
        output_by_target = dict()
        for project in pipeline.projects:
            for s1 in project.samples:
                for s2 in project.samples:
                    if s1 == s2:
                        continue
                    pair = Pair(s1, s2)
                    output_by_target[pair.unique_id] = \
                        self._queue_jobs_with_checks(pair)
        return output_by_target


class ProjectStage(Stage, ABC):
    """
    Project-level stage
    """
    @abstractmethod
    def get_expected_output(
        self, 
        project: 'Project'
    ) -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, project: 'Project', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> Dict[str, StageOutput]:
        output_by_target = dict()
        for project in pipeline.projects:
            output_by_target[project.unique_id] = \
                self._queue_jobs_with_checks(project)
        return output_by_target
        

class CohortStage(Stage, ABC):    
    """
    Entire cohort level stage
    """
    @abstractmethod
    def get_expected_output(
        self, 
        pipeline: 'Pipeline'
    ) -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, pipe: 'Pipeline', inputs: StageInput) -> StageOutput:
        pass

    def add_to_the_pipeline(self, pipeline: 'Pipeline') -> Dict[str, StageOutput]:
        return {
            pipeline.unique_id: 
            self._queue_jobs_with_checks(pipeline)
        }


# Global _pipeline pointer to make the `@stage` decorator work before 
# a Pipeline instance is defined (i.e. before `run_pipeline` is called)
_pipeline = None


def run_pipeline(*args, stages_in_order: List[StageDecorator], **kwargs):
    """
    This function should be called to trigger the pipeline logic (i.e. find samples,
    call Stages' add_jobs(), submit Batch). Sort of implements the singleton logic
    for the Pipeline class, so there can be no more than one pipeline at runtime
    """
    global _pipeline
    if _pipeline is None:
        _pipeline = Pipeline(*args, **kwargs)
    _pipeline.set_stages([cls(_pipeline) for cls in stages_in_order])
    _pipeline.submit_batch()


class Pipeline(Target):
    """
    Represents a Pipeline, and incapulates projects, samples and a Hail Batch object.

    Inherited from StageTarget for global (cohort-level) stage (of type CohortStage)
    """
    def __init__(
        self,
        analysis_project: str,
        name: str,
        output_version: str,
        namespace: Namespace,
        title: str,
        keep_scratch: bool = False,
        dry_run: bool = False,
        previous_batch_tsv_path: Optional[str] = None,
        previous_batch_id: Optional[str] = None,
        update_smdb_analyses: bool = False,
        check_smdb_seq_existence: bool = False,
        skip_samples_without_seq_input: bool = False,
        validate_smdb_analyses: bool = False,
        check_intermediate_existence: bool = True,
        hail_billing_project: Optional[str] = None,
        first_stage: Optional[str] = None,
        last_stage: Optional[str] = None,
        config: Optional[Dict] = None,
        input_projects: Optional[List[str]] = None,
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        force_samples: Optional[List[str]] = None,
        ped_files: Optional[List[str]] = None,
    ):
        super().__init__()
        self.analysis_project = Project(
            pipeline=self,
            name=analysis_project,
            namespace=namespace,
        )
        self.name = name
        self.output_version = output_version
        self.namespace = namespace
        self.check_intermediate_existence = check_intermediate_existence
        self.first_stage = first_stage
        self.last_stage = last_stage
        self.db = SMDB(
            self.analysis_project.name,
            do_update_analyses=update_smdb_analyses,
            do_check_seq_existence=check_smdb_seq_existence,
        )
        self.skip_samples_without_seq_input = skip_samples_without_seq_input
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
            f'https://{self.namespace}-web.populationgenomics.org.au/'
            f'{self.analysis_project.stack}/'
            f'{self.name}/'
            f'{self.output_version}'
        )
        self.local_tmp_dir = tempfile.mkdtemp()
        self.keep_scratch = keep_scratch
        self.dry_run = dry_run

        self.prev_batch_jobs = dict()
        if previous_batch_tsv_path is not None:
            assert previous_batch_id is not None
            self.prev_batch_jobs = PrevJob.parse(
                previous_batch_tsv_path,
                previous_batch_id,
                get_hail_bucket(self.tmp_bucket, keep_scratch),
            )

        self.config = config or {}

        self.b = setup_batch(
            title, 
            self.keep_scratch, 
            self.tmp_bucket, 
            self.analysis_project.name, 
            hail_billing_project
        )

        self._stages_dict: Dict[str, Stage] = dict()
        self.projects: List[Project] = []
        if input_projects:
            self._populate(
                input_projects=input_projects,
                namespace=self.namespace,
                skip_samples=skip_samples,
                only_samples=only_samples,
                ped_files=ped_files,
                forced_samples=force_samples,
            )

    def get_all_samples(self) -> List[Sample]:
        all_samples = []
        for proj in self.projects:
            all_samples.extend(proj.samples)
        return all_samples

    def get_all_sample_ids(self) -> List[str]:
        return [s.id for s in self.get_all_samples()]

    def submit_batch(self) -> None:
        if self.b:
            logger.info(f'Will submit {self.b.total_job_num} jobs:')
            for label, stat in self.b.labelled_jobs.items():
                logger.info(f'  {label}: {stat["job_n"]} for '
                            f'{len(stat["samples"])} samples')
            logger.info(f'  Other jobs: {self.b.other_job_num}')

            self.b.run(
                dry_run=self.dry_run,
                delete_scratch_on_exit=not self.keep_scratch,
                wait=False,
            )
        shutil.rmtree(self.local_tmp_dir)

    def _populate_projects(
        self,
        input_projects: List[str],
        namespace: Namespace,
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        forced_samples: Optional[List[str]] = None,
    ):
        samples_by_project = self.db.get_samples_by_project(
            project_names=input_projects,
            namespace=namespace,
            skip_samples=skip_samples,
            only_samples=only_samples,
        )
        for proj_name, sample_datas in samples_by_project.items():
            project = Project(
                pipeline=self,
                name=proj_name,
                namespace=namespace,
            )
            for s_data in sample_datas:
                s = Sample(
                    id=s_data['id'], 
                    external_id=s_data['external_id'],
                    project=project,
                )
                if forced_samples and s.id in forced_samples:
                    logger.info(f'Force rerunning sample {s.id} even if outputs exist')
                    s.forced = True
                project.samples.append(s)
            self.projects.append(project)

    def _populate_analysis(self):
        all_sample_ids = self.get_all_sample_ids()

        jc_analysis = self.db.find_joint_calling_analysis(
            sample_ids=all_sample_ids,
        )

        for project in self.projects:
            sample_ids = [s.id for s in project.samples]
            
            seqs_by_sid = self.db.find_seq_by_sid(sample_ids)

            cram_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='cram',
            )
            gvcf_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='gvcf',
            )

            for s in project.samples:
                s.seq_info = seqs_by_sid.get(s.id)
                s.analysis_by_type[AnalysisType.CRAM] = cram_per_sid.get(s.id)
                s.analysis_by_type[AnalysisType.GVCF] = gvcf_per_sid.get(s.id)

        self.analysis_by_type[AnalysisType.JOINT_CALLING] = jc_analysis

    def _populate_pedigree(self, ped_files: List[str]):
        sample_by_extid = dict()
        for s in self.get_all_samples():
            sample_by_extid[s.external_id] = s

        for i, ped_file in enumerate(ped_files):
            local_ped_file = join(self.local_tmp_dir, f'ped_file_{i}.ped')
            utils.gsutil_cp(ped_file, local_ped_file)
            with open(local_ped_file) as f:
                for line in f:
                    if 'Family.ID' in line:
                        continue
                    fields = line.strip().split('\t')[:6]
                    fam_id, sam_id, pat_id, mat_id, sex, phenotype = fields
                    if sam_id in sample_by_extid:
                        s = sample_by_extid[sam_id]
                        s.pedigree = PedigreeInfo(
                            fam_id=fam_id,
                            dad=sample_by_extid.get(pat_id),
                            mom=sample_by_extid.get(mat_id),
                            sex={
                                '1': Sex.MALE, 
                                '2': Sex.FEMALE,
                                'M': Sex.MALE,
                                'F': Sex.FEMALE,
                            }.get(sex, Sex.UNKNOWN),
                            phenotype=phenotype or '0',
                        )
        for project in self.projects:
            samples_with_ped = [s for s in project.samples if s.pedigree]
            logger.info(
                f'{project.name}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(project.samples)}'
            )

    def _populate(
        self,
        input_projects: List[str],
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        ped_files: Optional[List[str]] = None,
        namespace: Optional[Namespace] = None,
        forced_samples: List[str] = None
    ) -> None:
        """
        Finds input samples, analyses and sequences from the DB, 
        populates self.projects, adds pedigree information
        """
        self._populate_projects(
            input_projects=input_projects,
            skip_samples=skip_samples,
            only_samples=only_samples,
            namespace=namespace or self.namespace,
            forced_samples=forced_samples,
        )
        self._populate_analysis()
        if ped_files:
            self._populate_pedigree(ped_files)

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

    def set_stages(self, stages: List[Stage]):
        for stage in stages:
            if stage.name in self._stages_dict:
                raise ValueError(
                    f'Stage {stage.name} is already defined. Check your '
                    f'list for duplicates: {", ".join(s.name for s in stages)}'
                )
            self._stages_dict[stage.name] = stage

        first_stage_num, last_stage_num = self._validate_first_last_stage()

        # First round - checking which stages we require, even if they are skipped
        for i, (stage_name, stage) in enumerate(self._stages_dict.items()):
            if first_stage_num is not None and i < first_stage_num:
                stage.skipped = True
                stage.required = False
                logger.info(f'Skipping stage {stage_name}')

            for reqcls in stage.required_stages_classes:
                assert reqcls.__name__ in self._stages_dict, (
                    reqcls.__name__, list(self._stages_dict.keys())
                )
                reqstage = self._stages_dict[reqcls.__name__]
                reqstage.required = True
                if reqstage.skipped:
                    logger.info(
                        f'Stage {reqstage.name} is skipped, '
                        f'but the output will be required for the stage {stage.name}'
                    )
                stage.required_stages.append(reqstage)

            if last_stage_num and i > last_stage_num:
                stage.skipped = True
        
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
        return utils.can_reuse(fpath, overwrite=not self.check_intermediate_existence)
    
    @property
    def unique_id(self) -> str:
        return self.name
