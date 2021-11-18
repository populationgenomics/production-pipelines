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
    def queue_jobs(self, sample: Sample, inputs: StageResults) -> StageResults:
        expected_path = f'{sample.project.get_bucket()}/cram/{sample.id}.cram'
        job = align.bwa(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, paths=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.GVCF, requires_stages=CramStage)
class GvcfStage(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageResults) -> StageResults:
        cram_path = inputs.get_file_path(target=sample, stage=CramStage)
        expected_path = f'{sample.project.get_bucket()}/gvcf/{sample.id}.g.vcf.gz'
        job = haplotype_caller.produce_gvcf(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(sample, paths=expected_path, jobs=[job])

@stage(analysis_type=AnalysisType.JOINT_CALLING)
class JointCallingStage(CohortStage):
    def queue_jobs(self, pipeline: Pipeline, inputs: StageResults) -> StageResults:
        gvcf_by_sid = {
            s.id: inputs.get_file_path(stage=GvcfStage, target=s)
            for p in pipe.projects
            for s in p.samples
        }
        expected_path = ...
        job = make_joint_genotyping_jobs(b=self.pipe.b, ..., output_path=expected_path)
        return self.make_outputs(pipe, paths=expected_path, jobs=[job])

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
import sys
import tempfile
from dataclasses import dataclass, field
import logging
from enum import Enum
from os.path import join
from typing import List, Dict, Optional, Any, Tuple, Callable, cast, Type, Union
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
    Decorator to use with click when writing a pipeline script, e.g.:
    @click.command()
    @pipeline_click_options
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
            help='SM project name to write the intermediate/joint-calling analysis entries to',
        ),
        click.option(
            '--input-project',
            'input_projects',
            multiple=True,
            required=True,
            help='Only read samples that belong to the project(s). Can be set multiple times.',
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
            help='Before running a job, check for an intermediate output before submitting it, '
                 'and if it exists on a bucket, submit a [reuse] job instead. Works well with '
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
                 'the bucket. Set the analysis entry to "failure" if output doesn\'t exist'
        ),
        click.option(
            '--previous-batch-tsv',
            'previous_batch_tsv_path',
            help='A list of previous successful attempts from another batch, dumped from '
                 'from the Batch database (the "jobs" table joined on "job_attributes"). '
                 'If the intermediate output for a job exists in a previous attempt, '
                 'it will be passed forward, and a [reuse] job will be submitted.'
        ),
        click.option(
            '--previous-batch-id',
            'previous_batch_id',
            help='6-letter ID of the previous successful batch (corresponds to the directory '
                 'name in the batch logs. e.g. feb0e9 in '
                 'gs://cpg-seqr-main-tmp/hail/batch/feb0e9'
        )
    ]
    for opt in options:
        function = opt(function)
    return function


class StageResults:
    """
    Object of this class is passed to `add_jobs` method of a Stage, and can be used
    to query dependency files and jobs
    """
    def __init__(
        self, 
        pipe: 'Pipeline',
        sourcestage: Optional['Stage'] = None,
    ):
        self._pipe = pipe
        self._sourcestage = sourcestage
        self._fpath_dict_by_target_by_stagecls: Dict[str, Dict[str, Dict[str, str]]] = dict()
        self._resource_dict_by_target_by_stagecls: Dict[str, Dict[str, Dict[str, hb.Resource]]] = dict()
        self._jobs: List[Job] = []
        
    def append(self, other: 'StageResults'):
        """
        Merge in another StageDeps object
        """
        for stage, fpath_dict_by_target in other._fpath_dict_by_target_by_stagecls.items():
            for target, fpath_dict in fpath_dict_by_target.items():
                self.add_paths(fpath_dict, stage=stage, target=target)

        for stage, resource_by_target in other._resource_dict_by_target_by_stagecls.items():
            for target, resource in resource_by_target.items():
                self.add_resources(resource, stage=stage, target=target)

        self.add_jobs(other.get_jobs())

    def get_file_path_by_target(
        self, 
        stage: Optional[Callable] = None,
        id: Optional[str] = None,
    ) -> Dict[str, str]:
        if stage is not None:
            stage_name = stage.__name__
        else:
            assert self._sourcestage is not None
            stage_name = self._sourcestage.name
        return {
            trg: self._fpath_dict_by_target_by_stagecls[stage_name][trg][id]
            for trg in self._fpath_dict_by_target_by_stagecls[stage_name]
        }

    def get_file_path_dict_by_target(
        self, 
        stage: Optional[Callable] = None,
    ) -> Dict[str, Dict[str, str]]:
        return self._get(
            d=self._fpath_dict_by_target_by_stagecls,
            stage=stage,
        )

    def get_file_path_dict(
        self, 
        target: 'StageTarget',
        stage: Optional[Callable] = None,
    ) -> Dict[str, str]:
        return self.get_file_path(
            target=target,
            stage=stage,
        )

    def get_file_path(
        self, 
        target: 'StageTarget',
        stage: Optional[Callable] = None,
        id: Optional[str] = None,
    ) -> str:
        return self._get(
            d=self._fpath_dict_by_target_by_stagecls,
            target=target,
            stage=stage,
            id=id
        )

    def get_resource_by_target(
        self, 
        stage: Optional[Callable] = None,
        id: Optional[str] = None,
    ) -> Dict[str, hb.Resource]:
        if stage is not None:
            stage_name = stage.__name__
        else:
            assert self._sourcestage is not None
            stage_name = self._sourcestage.name
        return {
            trg: self._resource_dict_by_target_by_stagecls[stage_name][trg][id]
            for trg in self._fpath_dict_by_target_by_stagecls[stage_name]
        }

    def get_resource_dict_by_target(
        self, 
        stage: Optional[Callable] = None,
    ) -> Dict[str, Dict[str, hb.Resource]]:
        return self._get(
            d=self._resource_dict_by_target_by_stagecls,
            stage=stage,
        )

    def get_resource_dict(
        self, 
        stage: Optional[Callable] = None,
        target: Optional['StageTarget'] = None
    ) -> Dict[str, hb.Resource]:
        return self.get_resource(
            stage=stage,
            target=target
        )

    def get_resource(
        self, 
        stage: Optional[Callable] = None,
        target: Optional['StageTarget'] = None,
        id: Optional[str] = None,
    ) -> hb.Resource:
        if target is None:
            target = self._pipe
        self._get(
            d=self._resource_dict_by_target_by_stagecls,
            stage=stage,
            target=target,
            id=id,
        )
    
    def _get(
        self,
        d: Dict,
        target: Optional['StageTarget'] = None,
        stage: Optional[Callable] = None,
        id: Optional[str] = None,
    ):
        if stage is not None:
            stage_name = stage.__name__
        else:
            assert self._sourcestage is not None
            stage_name = self._sourcestage.name
        d = d[stage_name]
        if target:
            d = d[target.get_target_name()]
        return d[id]

    def get_jobs(self) -> List[Job]:
        return self._jobs

    def add_jobs(self, jobs: Optional[List[Job]]):
        if jobs:
            self._jobs.extend(jobs)

    def add_paths(
        self, 
        data: Union[str, Dict[str, str]],
        stage: Optional[str] = None,
        target: Optional[Union['StageTarget', str]] = None,
    ):
        """
        Add dependency file paths
        """
        self._add_file_data(data, stage=stage, is_resource=False, target=target)

    def add_resources(
        self, 
        data: Union[hb.Resource, Dict[str, hb.Resource]],
        stage: Optional[str] = None,
        target: Optional[Union['StageTarget', str]] = None,
    ):
        """
        Add dependency Hail Batch resources
        """
        self._add_file_data(data, stage=stage, is_resource=True, target=target)

    def _add_file_data(
        self,
        data: Union[Union[str, hb.Resource], Dict[str, Union[str, hb.Resource]]],
        is_resource: bool,
        stage: Optional[str] = None,
        target: Optional[Union['StageTarget', str]] = None,
    ):
        """
        Add dependency file paths or Hail Batch output resoruces.
        Can be one object or a dictionary.
        """
        if stage is None:
            assert self._sourcestage is not None
            stage = self._sourcestage.name
        if target is None:
            target = self._pipe
        if not isinstance(data, dict):
            data = {None: data}
        if is_resource:
            d = self._resource_dict_by_target_by_stagecls
        else:
            d = self._fpath_dict_by_target_by_stagecls
        if stage not in d:
            d[stage] = dict()
        if not isinstance(target, str):
            target = cast(StageTarget, target).get_target_name()
        d[stage][target] = data


# class StageLevel(Enum):
#     SAMPLE = 0
#     PROJECT = 0
#     COHORT = 0


def stage(
    _cls=None, 
    *,
    analysis_type: Optional[AnalysisType] = None, 
    requires_stages: Optional[Union[List[Type['Stage']], Type]] = None,
    # level: StageLevel = StageLevel.COHORT,
):
    """
    Implements a standard class decorator pattern with an optional argument.
    The goal is to make subclasses of Stage singleton classes, that would allow 
    only one object per subclass.
    """
    def decorator_stage(cls):
        @functools.wraps(cls)
        def wrapper_stage(pipe: 'Pipeline' = None):
            return cls(
                name=cls.__name__,
                pipe=pipe,
                analysis_type=analysis_type,
                requires_stages=requires_stages,
                # level=level,
            )
        return wrapper_stage

    if _cls is None:
        return decorator_stage
    else:
        return decorator_stage(_cls)


# def sample_stage(*args, **kwargs):
#     return stage(*args, **kwargs, level=StageLevel.SAMPLE)
# 
# 
# def project_stage(*args, **kwargs):
#     return stage(*args, **kwargs, level=StageLevel.PROJECT)
# 
# 
# def cohort_stage(*args, **kwargs):
#     return stage(*args, **kwargs, level=StageLevel.COHORT)


# def sample_stage(
#     _fun=None, 
#     *, 
#     analysis_type: Optional[AnalysisType] = None, 
#     requires_stages: Optional[Union[List[Type['Stage']], Type]] = None
# ):
#     """
#     Implements a standard Class decorator pattern with an optional argument.
#     The goal is to make subclasses of Stage singleton classes, that would allow 
#     only one object per subclass.
#     """
#     stage.instances = dict()  # stores instance per subclass
# 
#     def decorator_stage(cls):
#         @functools.wraps(cls)
#         def wrapper_stage(*args, **kwargs):
#             if cls not in stage.instances:
#                 stage.instances[cls] = cls(
#                     analysis_type=analysis_type,
#                     requires_stages=requires_stages
#                 )
#             return stage.instances[cls]
#         return wrapper_stage
# 
#     if _cls is None:
#         return decorator_stage
#     else:
#         return decorator_stage(_cls)


class Stage(ABC):
    """
    Abstract class for a pipeline stage
    """
    def __init__(
        self,
        pipe: 'Pipeline',
        name: str,
        requires_stages: Optional[Union[List[Callable], Type]] = None,
        analysis_type: Optional[AnalysisType] = None,
    ):
        self.name = name
        self.pipe: 'Pipeline' = pipe
        self.required_stages_classes: List[Callable] = []
        if requires_stages:
            if isinstance(requires_stages, list):
                self.required_stages_classes.extend(requires_stages)
            elif isinstance(requires_stages, Callable):
                self.required_stages_classes.append(requires_stages)
            else:
                logger.critical(
                    'requires_stages must be a Type or a list of Type objects that ' +
                    'are subclasses of Stage (SampleStage, ProjectStage, CohortStage)'
                )
                sys.exit(1)

        # Populated in pipeline.run(), after we know all stages
        self.required_stages: List[Stage] = []
        
        # If analysis type is defined, it will be used to find and reuse existing
        # outputs from the SMDB
        self.analysis_type = analysis_type

        # Populated when the stage is added by calling `add_to_the_pipeline`
        self.results: StageResults = StageResults(pipe=self.pipe, sourcestage=self)

        # self.active=True means that the stage wasn't skipped and jobs were added
        # into the pipeline, and we shoule expect target.ouptut_by_stage 
        # to be populated. Otherwise, self.get_expected_output() should work.
        self.skipped = False  
        # self.required=True means that is required for another active stage,
        # even if it was skipped
        self.required = True

    def get_name(self):
        return self.name

    @staticmethod
    @abstractmethod
    def get_expected_output(target: 'StageTarget') -> Optional[Union[str, Dict[str, str]]]:
        """
        Path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """
        pass

    @abstractmethod
    def queue_jobs(
        self,
        target: 'StageTarget',
        inputs: StageResults,
    ) -> StageResults:
        """
        Implements logic of the Stage: creates Batch jobs that do the processing.
        Assumes that all the household work is done: checking missing inputs
        from requried stages, checking the SMDB, checking for possible reuse of 
        existing outputs.
        """
        pass

    @abstractmethod
    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call `output = add_jobs(target, input)` for each target in the pipeline
        """
        
    def _add_to_the_pipeline_for_target(self, target: 'StageTarget') -> StageResults:
        if not self.skipped:
            logger.info(f'{self.get_name()}: adding jobs for {target.get_target_name()}')
            return self._add_or_reuse_jobs(target)
        elif self.required:
            reuse_paths = self._check_if_can_reuse(target)
            if not reuse_paths:
                logger.critical(
                    f'Stage {self.name} is required, but skipped and '
                    f'cannot reuse outputs for {target.get_target_name()}'
                )
                sys.exit(1)
            return self.make_outputs(target, paths=reuse_paths) 
        else:
            # Stager is not needed, returning empty outputs
            return self.make_outputs(target)

    # def append_to_inputs_from_prev_stage(
    #     self,
    #     inputs: StageResults,
    #     prev_stage: 'Stage',
    #     prev_target: 'StageTarget',
    # ):
    #     """
    #     Populates output genereated by `prev_stage` into `inputs` for each 
    #     `prev_target`.
    #     
    #     Sort of a multiple dispatch. 
    #     
    #     Collecting dependencies paths. Will return a path when the dependency is the
    #     same or a higher-level stage (e.g. Project-level to a Sample-level), or will
    #     return a dictionary when the dependency is of a more fine-grained level:
    #     (e.g. outputs of a Sample-level stage for a Project-level stage)
    #     """
    #     if (
    #         issubclass(self.__class__, SampleStage) and issubclass(prev_stage.__class__, SampleStage) or
    #         issubclass(self.__class__, ProjectStage) and issubclass(prev_stage.__class__, ProjectStage) or
    #         issubclass(self.__class__, CohortStage) and issubclass(prev_stage.__class__, CohortStage)
    #     ):
    #         inputs.add_paths(prev_stage.find_output(prev_target))
    #         inputs.add_jobs(prev_target.jobs_by_stagecls.get(prev_stage.__class__))
    # 
    #     elif issubclass(self.__class__, ProjectStage) and issubclass(prev_stage.__class__, SampleStage):
    #         dep_path_by_id = dict()
    #         dep_jobs = []
    #         for s in cast(Project, prev_target).samples:
    #             dep_path_by_id[cast(StageTarget, s)] = prev_stage.find_output(s)
    #             dep_jobs.extend(s.jobs_by_stagecls.get(prev_stage.__class__, []))
    #         return dep_path_by_id, dep_jobs
    # 
    #     elif issubclass(self.__class__, CohortStage) and issubclass(prev_stage.__class__, SampleStage):
    #         dep_path_by_id = dict()
    #         dep_jobs = []
    #         for p in cast(Pipeline, prev_target).projects:
    #             for s in p.samples:
    #                 dep_path_by_id[cast(StageTarget, s)] = prev_stage.find_output(s)
    #                 dep_jobs.extend(s.jobs_by_stagecls.get(prev_stage.__class__, []))
    #         return dep_path_by_id, dep_jobs
    # 
    #     elif issubclass(self.__class__, CohortStage) and issubclass(prev_stage.__class__, ProjectStage):
    #         dep_path_by_id = dict()
    #         dep_jobs = []
    #         for p in cast(Pipeline, prev_target).projects:
    #             dep_path_by_id[cast(StageTarget, p)] = prev_stage.find_output(p)
    #             dep_jobs.extend(p.jobs_by_stagecls.get(prev_stage.__class__, []))
    #         return dep_path_by_id, dep_jobs
    # 
    #     elif issubclass(self.__class__, SampleStage) and issubclass(prev_stage.__class__, ProjectStage):
    #         sample = cast(Sample, prev_target)
    #         return (
    #             {self.pipe: prev_stage.find_output(sample.project)},
    #             sample.project.jobs_by_stagecls.get(prev_stage.__class__, [])
    #         )
    # 
    #     elif issubclass(prev_stage.__class__, CohortStage):
    #         return (
    #             {self.pipe: prev_stage.find_output(self.pipe)},
    #             self.pipe.jobs_by_stagecls.get(prev_stage.__class__, [])
    #         )
    #     else:
    #         logger.critical(f'Unknown stage types: {type(self), type(prev_stage)}')
    #         sys.exit(1)

    def _add_jobs(self, target: 'StageTarget') -> StageResults:
        """
        Checks the existence of the depenencies, and calls 
        the public `outputs = self.queue_jobs(target, inputs)`
        """
        # Building the inputs for the job here
        inputs = StageResults(pipe=self.pipe)

        for prev_stage in self.required_stages:
            # d = prev_stage.results.get_file_path_dict_by_target()
            # for target, file_path_dict in d.items():
            #     for k, file_path in file_path_dict.items():
            #         if not self.pipe.check_intermediate_existence:
            # 
            #         if not k:
            
            inputs.append(prev_stage.results)

        return self.queue_jobs(target, inputs)
            
            # self.append_to_inputs_from_prev_stage(
            #     inputs=inputs,
            #     prev_stage=prev_stage, 
            # )

            # target_name = target.get_target_name()
            # target_name = f' (on {target_name})' if target_name else ''
            # 
            # if isinstance(list(deppath_by_deptarget.keys())[0], Pipeline):
            #     # When we deal with stages of the same level 
            #     # (Stage-Stage, Project-Project, Cohort-Cohort), they will be indexed
            #     # by a Pipeline object
            #     deppathdata = list(deppath_by_deptarget.values())[0]
            #     if deppathdata is None:
            #         logger.critical(
            #             f'Cannot find required output(s) '
            #             f'for stage {self.get_name()}{target_name} '
            #             f'from stage {prev_stage.get_name()}. Exiting'
            #         )
            #         sys.exit(1)
            #     else:
            #         inputs.add_paths(deppathdata, stagecls=reqcls)
            # else:
            #     targets_with_not_found_deps: List['StageTarget'] = []
            #     for deptarget, deppathdata in deppath_by_deptarget.items():
            #         if deppathdata is None:
            #             targets_with_not_found_deps.append(deptarget)
            #         else:
            #             inputs.add_paths(deppathdata, target=deptarget, stagecls=reqcls)
            # 
            #     if targets_with_not_found_deps:
            #         logger.critical(
            #             f'Cannot find {len(targets_with_not_found_deps)} required outputs '
            #             f'for stage {self.get_name()}{target_name} '
            #             f'from stage {prev_stage.get_name()} on '
            #             f'{"|".join(str(t.get_target_name()) for t in targets_with_not_found_deps)}. '
            #             f'Exiting'
            #         )
            #         sys.exit(1)
            # 
            # inputs.add_jobs(dep_jobs, stagecls=reqcls)

        # stage_results: StageResults = self.queue_jobs(target, inputs)
        # res_file_data = stage_results.get_file_path_dict()
        # res_resource_data = stage_results.get_resource_dict()
        # res_jobs = stage_results.get_jobs()
        # if res_jobs:
        #     target.jobs_by_stagecls[self.__class__] = res_jobs
        # return stage_results
    
    def _check_if_can_reuse(self, target: 'StageTarget') -> Optional[Union[Dict[str, str], str]]:
        expected_paths = self.get_expected_output(target)

        if self.analysis_type is not None:
            if isinstance(expected_paths, dict):
                logger.error(
                    f'get_expected_output() returns a dict, won\'t check the SMDB for '
                    f'{self.name} on {target.get_target_name()}'
                )
            expected_path = cast(str, expected_paths)
            found_path = self._check_smdb_analysis(target, expected_path)
            if found_path:
                return found_path
        if expected_paths and self.pipe.validate_smdb_analyses:
            if utils.file_exists(expected_paths.values()):
                return expected_paths
        return None

    def _reuse_jobs(
        self, 
        target: 'StageTarget', 
        found_paths: Union[str, Dict[str, str]]
    ) -> StageResults:
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
            paths=found_paths,
            jobs=[self.pipe.b.new_job(f'{self.get_name()} [reuse]', attributes)]
        )

    def _add_or_reuse_jobs(self, target: 'StageTarget') -> StageResults:
        reuse_paths = self._check_if_can_reuse(target)
        if reuse_paths:
            return self._reuse_jobs(target, reuse_paths)
        else:
            return self._add_jobs(target)

    def _check_smdb_analysis(
        self, 
        target: 'StageTarget',
        expected_path: str,
    ) -> Optional[str]:
        """
        Check if SMDB already has analysis, and invalidate it if the
        output doesn't exist
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
    
    def make_outputs(
        self, 
        target: 'StageTarget',
        resources: Optional[Union[Dict[str, hb.Resource], hb.Resource]] = None, 
        paths: Optional[Union[str, Dict[str, str]]] = None,
        jobs: Optional[List[Job]] = None
    ) -> StageResults:
        """
        Builds a StageDeps object to return from a stage's add_jobs()
        """
        deps = StageResults(pipe=self.pipe, sourcestage=self)
        if jobs:
            deps.add_jobs(jobs)
        if resources:
            deps.add_resources(resources, target=target)
        if paths:
            deps.add_paths(paths, target=target) 
        return deps
        

class SampleStage(Stage, ABC):
    """
    Sample-level stage. run() takes sample and project name as input
    """
    @abstractmethod
    def queue_jobs(
        self, 
        target: 'Sample',
        inputs: StageResults,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass

    def _add_or_reuse_jobs(self, sample: 'Sample') -> StageResults:
        if sample.id in self.pipe.force_samples:
            return self._add_jobs(sample)
        else:
            return super()._add_or_reuse_jobs(sample)        

    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call `output = add_jobs(target, input)` for each target in the pipeline
        """
        for project in pipe.projects:
            for sample in project.samples:
                sample_res = self._add_to_the_pipeline_for_target(sample)
                self.results.append(sample_res)


class ProjectStage(Stage, ABC):
    """
    Project-level stage. run() takes project name as input
    """
    @abstractmethod
    def queue_jobs(
        self,
        target: 'Project',
        inputs: StageResults,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass

    def add_to_the_pipeline(self, pipe: 'Pipeline') -> StageResults:
        """
        Call `output = add_jobs(target, input)` for each target in the pipeline
        """
        for project in pipe.projects:
            self.results.append(self._add_to_the_pipeline_for_target(project))


class CohortStage(Stage, ABC):    
    """
    Entire cohort level stage
    """
    @abstractmethod
    def queue_jobs(
        self, 
        target: 'Pipeline',
        inputs: StageResults,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass

    def add_to_the_pipeline(self, pipe: 'Pipeline') -> StageResults:
        """
        Call `output = add_jobs(target, input)` for each target in the pipeline
        """
        self.results = self._add_to_the_pipeline_for_target(pipe)


class StageTarget:
    """
    Sample, Project, Cohort inherit from this class
    """
    def __init__(self):
        # From SMDB Analysis entries:
        self.analysis_by_type: Dict[AnalysisType, Analysis] = dict()
    
    @abstractmethod
    def get_target_name(self) -> Optional[str]:
        pass


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


@dataclass(init=False)
class Sample(StageTarget):
    """
    Corresponds to one Sample entry in the SMDB
    """
    id: str
    external_id: str
    project: 'Project' = field(repr=lambda p: p.name)
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

    def get_target_name(self) -> Optional[str]:
        return self.id


@dataclass(init=False)
class Project(StageTarget):
    """
    Represents a CPG stack (aka dataset) or a pair of main-test projects in the SDMB) 
    """
    stack: str  # in stack terms: e.g. can be seqr but _not_ seqr-test
    name: str   # in smdb terms: e.g. can be seqr, seqr-test
    samples: List[Sample]
    pipeline: 'Pipeline'
    is_test: bool

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

    def __init__(
        self, 
        pipeline: 'Pipeline',
        name: str, 
        namespace: Namespace
    ):
        super().__init__()
        self.pipeline = pipeline
        self.samples = []
        self.is_test = namespace != Namespace.MAIN
        if name.endswith('-test'):
            self.is_test = True
            self.stack = name[:-len('-test')]
        else:
            self.stack = name
        
        self.name = self.stack
        if self.is_test:
            self.name = self.stack + '-test'
            
    def get_target_name(self) -> Optional[str]:
        return self.name


_pipeline: 'Pipeline' = None


def run_pipeline(*args, stages_in_order: List[Callable], **kwargs):
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


class Pipeline(StageTarget):
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
        config: Any = None,
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
        self.force_samples = force_samples or []
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
        self.prev_batch_jobs = PrevJob.parse(
            previous_batch_tsv_path,
            previous_batch_id,
            get_hail_bucket(self.tmp_bucket, keep_scratch),
        ) if previous_batch_tsv_path else dict()
        self.config = config

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
                logger.info(f'  {label}: {stat["job_n"]} for {len(stat["samples"])} samples')
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
                project.samples.append(Sample(
                    id=s_data['id'], 
                    external_id=s_data['external_id'],
                    project=project,
                ))
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
            if stage.get_name() in self._stages_dict:
                logger.critical(
                    f'Stage {stage.get_name()} is already defined. Check your '
                    f'list for duplicates: {", ".join(s.get_name() for s in stages)}'
                )
                sys.exit(1)
            self._stages_dict[stage.get_name()] = stage

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
                        f'Stage {reqstage.get_name()} is skipped, '
                        f'but the output will be required for the stage {stage.get_name()}'
                    )
                stage.required_stages.append(reqstage)

            logger.info(f'')
            if last_stage_num and i > last_stage_num:
                stage.skipped = True

        # Second round - actually adding jobs from the stages
        for i, (stagecls, stage) in enumerate(self._stages_dict.items()):
            if not stage.skipped:
                logger.info(f'*' * 60)
                logger.info(f'Stage {stage.get_name()}')

            if stage.required:
                stage.add_to_the_pipeline(self)

            if not stage.skipped:
                logger.info(f'')
                if last_stage_num and i >= last_stage_num:
                    logger.info(f'Last stage is {stage.get_name()}, stopping here')
                    break

    def can_reuse(self, fpath: Optional[str]) -> bool:
        """
        Checks if the fpath exists, 
        but always returns False if not check_intermediate_existence
        """
        return utils.can_reuse(fpath, overwrite=not self.check_intermediate_existence)

    def get_target_name(self) -> Optional[str]:
        return None
