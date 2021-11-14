import os
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
import logging
from enum import Enum
from os.path import join
from typing import List, Dict, Optional, Any, Union, Tuple
from abc import ABC, abstractmethod

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import utils
from cpg_pipes.utils import Namespace, AnalysisType
from cpg_pipes.smdb import SMDB, Analysis, Sequence
from cpg_pipes.hailbatch import AlignmentInput, PrevJob

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Batch(hb.Batch):
    """
    Batch wrapper that registers job statistics
    """
    def __init__(self, name, backend, *args, **kwargs):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry
        self.labelled_jobs = dict()
        self.other_job_num = 0
        self.total_job_num = 0

    def new_job(
        self,
        name: Optional[str] = None,
        attributes: Optional[Dict[str, str]] = None,
        **kwargs,
    ) -> Job:
        """
        Adds job to the Batch, and also registers it in `self.job_stats` for
        statistics.
        """
        if not name:
            logger.critical('Error: job name must be defined')
        
        attributes = attributes or dict()
        project = attributes.get('project')
        sample = attributes.get('sample')
        samples = attributes.get('samples')
        label = attributes.get('label', name)

        if sample and project:
            name = f'{project}/{sample}: {name}'
        elif project:
            name = f'{project}: {name}'

        if label and (sample or samples):
            if label not in self.labelled_jobs:
                self.labelled_jobs[label] = {'job_n': 0, 'samples': set()}
            self.labelled_jobs[label]['job_n'] += 1
            self.labelled_jobs[label]['samples'] |= (samples or {sample})
        else:
            self.other_job_num += 1
        self.total_job_num += 1
        j = super().new_job(name, attributes=attributes)
        return j
    

class PipelineStage(ABC):
    def __init__(
        self, 
        pipe: 'Pipeline', 
        analysis_type: AnalysisType = AnalysisType.CUSTOM,
        requires_stages: Optional[List] = None
    ):
        """
        requires_stages: list of stage names
        """
        self.pipe = pipe
        self.active = False  
        # self.active=True means taht the stage wasn't skipped and jobs were added
        # into the pipeline, and we shoule expect target.ouptut_by_stage 
        # to be populated. Otherwise, self.get_expected_output() should work.
        self.analysis_type = analysis_type
        self.requires_stages = requires_stages or []

    @classmethod
    def get_name(cls):
        return cls.__name__.replace('Stage', '')

    @abstractmethod
    def get_expected_output(self, *args):
        pass

    @abstractmethod
    def add_jobs(
        self,
        target: Union['Sample', 'Project', 'Pipeline'],
        dep_paths_by_stage: Optional[Dict[str, Union[Dict[str, str], str]]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass
    
    @abstractmethod
    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call add_jobs(target) for each target in the pipeline
        """
        pass
    
    def get_deps(
        self, 
        dep: 'PipelineStage',  # class inherited from PipelineStage
        target: Union['Sample', 'Project', 'Pipeline'],
    ) -> Tuple[Union[Dict[str, str], str], List[Job]]:
        """
        Sort of a multiple dispatch. 
        
        Collecting dependencies paths. Will return a path when the dependency is the
        same or a higher-level stage (e.g. Project-level to a Sample-level), or will
        return a dictionary when the dependency is of a more fine-grained level:
        (e.g. outputs of a Sample-level stage for a Project-level stage)
        """
        if (
            issubclass(self.__class__, SampleStage) and issubclass(dep.__class__, SampleStage) or
            issubclass(self.__class__, ProjectStage) and issubclass(dep.__class__, ProjectStage) or
            issubclass(self.__class__, CohortStage) and issubclass(dep.__class__, CohortStage)
        ):
            return(
                dep.find_output(target),
                target.jobs_by_stage.get(dep.get_name(), [])
            )

        elif issubclass(self.__class__, ProjectStage) and issubclass(dep.__class__, SampleStage):
            dep_path_by_id = dict()
            dep_jobs = []
            project: Project = target
            for s in project.samples:
                dep_path_by_id[s.id] = dep.find_output(s)
                dep_jobs.extend(s.jobs_by_stage.get(dep.get_name(), []))
            return dep_path_by_id, dep_jobs

        elif issubclass(self.__class__, CohortStage) and issubclass(dep.__class__, SampleStage):
            dep_path_by_id = dict()
            dep_jobs = []
            pipe: 'Pipeline' = target
            for p in pipe.projects:
                for s in p.samples:
                    dep_path_by_id[s.id] = dep.find_output(s)
                    dep_jobs.extend(s.jobs_by_stage.get(dep.get_name(), []))
            return dep_path_by_id, dep_jobs

        elif issubclass(self.__class__, CohortStage) and issubclass(dep.__class__, ProjectStage):
            dep_path_by_id = dict()
            dep_jobs = []
            pipe: 'Pipeline' = target
            for p in pipe.projects:
                dep_path_by_id[p.name] = dep.find_output(p)
                dep_jobs.extend(p.jobs_by_stage.get(dep.get_name(), []))
            return dep_path_by_id, dep_jobs

        elif issubclass(self.__class__, SampleStage) and issubclass(dep.__class__, ProjectStage):
            sample: Sample = target
            return (
                dep.find_output(sample.project),
                sample.project.jobs_by_stage.get(dep.get_name(), [])
            )

        elif issubclass(dep.__class__, CohortStage):
            return (
                dep.find_output(self.pipe),
                self.pipe.jobs_by_stage.get(dep.get_name(), [])
            )
        else:
            logger.critical(f'Unknown stage types: {type(self), type(dep)}')
            sys.exit(1)

    def _add_jobs(
        self,
        target: Union['Sample', 'Project', 'Pipeline'],
    ):
        dep_paths_by_stage = dict()
        all_dep_jobs = []
        for depcls in self.requires_stages:
            if depcls.get_name() not in self.pipe.stage_by_name:
                logger.critical(f'Stage {depcls.get_name()} is not found')
            dep = self.pipe.stage_by_name.get(depcls.get_name())
            dep_paths, dep_jobs = self.get_deps(dep, target)
            all_dep_jobs.extend(dep_jobs)
            target_name = self.target_name(target)
            target_name = f"/{target_name}" if target_name else ""
            if isinstance(dep_paths, dict):
                not_found_deps = dict()
                for id, dep_path in dep_paths.items():
                    if not dep_path:
                        not_found_deps[id] = dep_path
                if not_found_deps:
                    logger.critical(
                        f'Cannot find {len(not_found_deps)} required outputs '
                        f'for stage {self.get_name()}{target_name} '
                        f'from stage {dep.get_name()}/'
                        f'{"|".join(not_found_deps.keys())}. Exiting'
                    )
                    sys.exit(1)
            else:
                if not dep_paths:
                    target_name = self.target_name(target)
                    logger.critical(
                        f'Cannot find required output '
                        f'for stage {self.get_name()}/{target_name} '
                        f'from stage {dep.get_name()}. Exiting'
                    )
                    sys.exit(1)
            dep_paths_by_stage[dep.get_name()] = dep_paths

        res_path, res_jobs = self.add_jobs(
            target,
            dep_paths_by_stage=dep_paths_by_stage,
            dep_jobs=all_dep_jobs,
        )
        if res_path:
            target.output_by_stage[self.get_name()] = res_path
        if res_jobs:
            target.jobs_by_stage[self.get_name()] = res_jobs

    def _reuse_jobs(
        self,
        target: Union['Sample', 'Project', 'Pipeline'],
        found_path: str,
    ):
        target.output_by_stage[self.get_name()] = found_path
        attributes = {}
        if isinstance(target, Sample):
            attributes['sample'] = target.id
            attributes['project'] = target.project.name
        if isinstance(target, Project):
            attributes['sample'] = target.name
        target.jobs_by_stage[self.get_name()] = [
            self.pipe.b.new_job(f'{self.get_name()} [reuse]', attributes)
        ]
        
    def _add_or_reuse_jobs(
        self, 
        target: Union['Sample', 'Project', 'Pipeline'],
    ):
        expected_path = self.get_expected_output(target)
        found_path = self._check_smdb_analysis(target, expected_path)
        if found_path:
            self._reuse_jobs(target, found_path)
        elif expected_path and self.pipe.validate_smdb_analyses and utils.file_exists(expected_path):
            self._reuse_jobs(target, expected_path)
        else:
            self._add_jobs(target)
    
    def find_output(
        self,
        target: Union['Sample', 'Project', 'Pipeline'],
    ) -> Optional[str]:
        """
        Called by the stage that depends on this one. This stage can be either active 
        (jobs were added), or skipped (in this case we can only return expected output)
        """
        if self.active:
            return target.output_by_stage.get(self.get_name())
        else:
            return self.get_expected_output(target)

    # def find_output(
    #     self,
    #     target: Union['Sample', 'Project', 'Pipeline'],
    #     expected_path: str,
    #     only_check_smdb: bool = False,
    #     only_check_bucket: bool = False,
    #     validate_smdb: bool = False,
    # ) -> Optional[str]:
    #     if not only_check_bucket:
    #         analysis = target.analysis_by_type.get(self.analysis_type)
    #         if only_check_smdb and analysis is None:
    #             return None
    # 
    #     if self.pipe.validate_smdb_analyses:
    #         if isinstance(target, Sample):
    #             sample: Sample = target
    #             sample_ids = [sample.id]
    #         elif isinstance(target, Project):
    #             project: Project = target
    #             sample_ids = [s.id for s in project.samples]
    #         else:
    #             pipe: Pipeline = target
    #             sample_ids = pipe.get_all_sample_ids()
    # 
    #         found_path = self.pipe.db.process_existing_analysis(
    #             sample_ids=sample_ids,
    #             completed_analysis=analysis,
    #             analysis_type=self.analysis_type,
    #             expected_output_fpath=expected_path,
    #         )
    #     else:
    #         found_path = analysis.output
    # 
    #     target.output_by_stage[self.get_name()] = found_path
    #     return found_path

    def _check_smdb_analysis(
        self, 
        target: Union['Sample', 'Project', 'Pipeline'],
        expected_path: str,
    ) -> Optional[str]:
        """
        Check if SMDB already has analysis, and invalidate it if the
        output doesn't exist
        """
        analysis = target.analysis_by_type.get(self.analysis_type)
        if not analysis:
            return None

        if self.pipe.validate_smdb_analyses:
            if isinstance(target, Sample):
                sample: Sample = target
                sample_ids = [sample.id]
            elif isinstance(target, Project):
                project: Project = target
                sample_ids = [s.id for s in project.samples]
            else:
                pipe: Pipeline = target
                sample_ids = pipe.get_all_sample_ids()

            found_path = self.pipe.db.process_existing_analysis(
                sample_ids=sample_ids,
                completed_analysis=analysis,
                analysis_type=self.analysis_type,
                expected_output_fpath=expected_path,
            )
        else:
            found_path = analysis.output

        target.output_by_stage[self.get_name()] = found_path
        return found_path
    
    @abstractmethod
    def target_name(self, target) -> Optional[str]:
        """
        Descriptive name of the stage target
        """
        pass
    
    
class SampleStage(PipelineStage, ABC):
    """
    Sample-level stage. run() takes sample and project name as input
    """
    @abstractmethod
    def add_jobs(
        self, 
        target: 'Sample',
        dep_paths_by_stage: Optional[Dict[str, Union[Dict[str, str], str]]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass
    
    def _add_or_reuse_jobs(self, sample: 'Sample'):
        if sample.id in self.pipe.force_samples:
            self._add_jobs(sample)
        else:
            super()._add_or_reuse_jobs(sample)        

    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call _add_jobs(target) for each target in the pipeline
        """
        for project in pipe.projects:
            logger.info(f'{self.get_name()}: adding jobs for project {project.name}')
            for sample in project.samples:
                logger.info(f'{self.get_name()}: adding jobs for {project.name}/{sample.id}')
                self._add_or_reuse_jobs(sample)     
    
    def target_name(self, sample: 'Sample') -> Optional[str]:
        """
        Descriptive name of the stage target
        """
        return f'{sample.project.name}/{sample.id}'


class ProjectStage(PipelineStage, ABC):
    """
    Project-level stage. run() takes project name as input
    """
    @abstractmethod
    def add_jobs(
        self,
        target: 'Project',
        dep_paths_by_stage: Optional[Dict[str, Union[Dict[str, str], str]]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass

    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call _add_jobs(target) for each target in the pipeline
        """
        for project in pipe.projects:
            logger.info(f'{self.get_name()}: adding jobs for project {project.name}')
            self._add_or_reuse_jobs(project)

    def target_name(self, project: 'Project') -> Optional[str]:
        """
        Descriptive name of the stage target
        """
        return f'{project.name}'
    

class CohortStage(PipelineStage, ABC):    
    """
    Entire cohort level stage
    """
    @abstractmethod
    def add_jobs(
        self, 
        target: 'Pipeline',
        dep_paths_by_stage: Optional[Dict[str, Union[Dict[str, str], str]]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        """
        Implements logic of the Stage: adds Batch jobs that do the processing.
        Assumes that all the household work is done - checking missing inputs,
        checking the SMDB, checking for possible reuse of existing outputs, etc
        """
        pass

    def add_to_the_pipeline(self, pipe: 'Pipeline'):
        """
        Call _add_jobs(target) for each target in the pipeline
        """
        logger.info(f'{self.get_name()}: adding jobs')
        self._add_or_reuse_jobs(pipe)

    def target_name(self, pipe: 'Pipeline') -> Optional[str]:
        """
        Descriptive name of the stage target
        """
        return None


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


@dataclass
class Sample:
    """
    Represents a Sample.
    """
    id: str
    external_id: str
    project: 'Project' = field(repr=lambda p: p.name)
    alignment_input: Optional[AlignmentInput] = None
    good: bool = True
    seq_info: Optional[Sequence] = None
    # From SMDB Analysis entries
    analysis_by_type: Dict[AnalysisType, Analysis] = field(default_factory=lambda: dict())
    # Other outputs, not supported by DB.
    # Project-level outputs are specified in Project class,
    # Cohort-level outputs are specified in Pipeline class
    output_by_stage: Dict[str, str] = field(default_factory=lambda: dict())
    jobs_by_stage: Dict[str, List[Job]] = field(default_factory=lambda: dict())

    pedigree: Optional[PedigreeInfo] = None

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


# TODO: potentially move get_sample_ids, jobs_by_stage, output_bu_stage
# into StageTarget, and pass StageTarget to stages instead of Sample or Project
# Not sure though what to do with analysis_by_type

@dataclass
class StageTarget(ABC):
    @abstractmethod
    def get_sample_ids(self):
        return ''


@dataclass
class SampleStageTarget(ABC):
    pass


@dataclass(init=False)
class Project:
    """
    Represents a CPG dataset/project.
    """
    stack: str  # in stack terms: can be seqr but not seqr-test
    name: str   # in smdb terms: can be seqr, seqr-test
    samples: List[Sample]
    pipeline: 'Pipeline'
    is_test: bool
    # From SMDB Analysis entries
    analysis_by_type: Dict[AnalysisType, Analysis]
    # Outputs not supported by DB
    # Sample-level outputs are specified and Sample class,
    # Cohort-level outputs are specified in Pipeline class
    output_by_stage: Dict[str, str]
    jobs_by_stage: Dict[str, Job]

    def get_bucket(self):
        return f'gs://cpg-{self.stack}-{self.pipeline.output_suf}'

    def get_tmp_bucket(self):
        return f'gs://cpg-{self.stack}-{self.pipeline.output_suf}-tmp'
    
    def __init__(
        self, 
        name: str, 
        pipeline: 'Pipeline',
        namespace: Namespace
    ):
        self.samples = []
        self.pipeline = pipeline
        self.is_test = namespace != Namespace.MAIN
        if name.endswith('-test'):
            self.is_test = True
            self.stack = name[:-len('-test')]
        else:
            self.stack = name
        
        self.name = self.stack
        if self.is_test:
            self.name = self.stack + '-test'
            
        self.analysis_by_type = dict()
        self.output_by_stage = dict()
        self.jobs_by_stage = dict()


class Pipeline:
    """
    Represents a processing pipeline, and incapsulates the Batch object, Batch jobs,
    samples and projects.
    """
    def __init__(
        self,
        analysis_project: str,
        name: str,
        output_version: str,
        namespace: Namespace,
        title: str,
        keep_scratch: bool = False,
        previous_batch_tsv_path: Optional[str] = None,
        previous_batch_id: Optional[str] = None,
        smdb_update_analyses: bool = False,
        smdb_check_seq_existence: bool = False,
        skip_samples_without_seq_input: bool = False,
        validate_smdb_analyses: bool = False,
        check_intermediate_existence: bool = True,
        hail_billing_project: Optional[str] = None,
        first_stage: Optional[str] = None,
        last_stage: Optional[str] = None,
        config: Any = None,
        input_projects: Optional[List[str]] = None,
        skip_samples: Optional[List[str]] = None,
        force_samples: Optional[List[str]] = None,
        ped_files: Optional[List[str]] = None,
    ):
        self.analysis_project = Project(
            name=analysis_project,
            pipeline=self,
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
            do_update_analyses=smdb_update_analyses,
            do_check_seq_existence=smdb_check_seq_existence,
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
            self.analysis_project, 
            hail_billing_project
        )

        # Found global SMDB analysis entries by type
        self.analysis_by_type: Dict[AnalysisType, Analysis] = dict()
        # Cohort-level outputs, not supported by SMDB.
        # Sample- and project-level outputs are specified in Sample and Project classes
        self.output_by_stage: Dict[str, str] = dict()
        self.jobs_by_stage: Dict[str, Job] = dict()

        self.stage_by_name = Dict

        self.projects: List[Project] = []
        if input_projects:
            self.populate(
                input_projects=input_projects,
                namespace=self.namespace,
                skip_samples=skip_samples,
                ped_files=ped_files,
            )

    def get_all_samples(self) -> List[Sample]:
        all_samples = []
        for proj in self.projects:
            all_samples.extend(proj.samples)
        return all_samples
    
    def get_all_sample_ids(self) -> List[str]:
        return [s.id for s in self.get_all_samples()]

    def run(self, dry_run: bool = False) -> None:
        if self.b:
            logger.info(f'Will submit {self.b.total_job_num} jobs:')
            for label, stat in self.b.labelled_jobs.items():
                logger.info(f'  {label}: {stat["job_n"]} for {len(stat["samples"])} samples')
            logger.info(f'  Other jobs: {self.b.other_job_num}')

            self.b.run(
                dry_run=dry_run,
                delete_scratch_on_exit=not self.keep_scratch,
                wait=False,
            )
        shutil.rmtree(self.local_tmp_dir)

    def _populate_projects(
        self,
        input_projects: List[str],
        namespace: Namespace,
        skip_samples: Optional[List[str]] = None,
    ):
        samples_by_project = self.db.get_samples_by_project(
            project_names=input_projects,
            namespace=namespace,
            skip_samples=skip_samples,
        )
        for proj_name, sample_datas in samples_by_project.items():
            project = Project(
                name=proj_name,
                pipeline=self,
                namespace=namespace,
            )
            for s_data in sample_datas:
                project.samples.append(Sample(
                    id=s_data['id'], 
                    external_id=s_data['external_id'],
                    project=project,
                    analysis_by_type=dict(),
                    output_by_stage=dict(),
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

    def populate(
        self,
        input_projects: List[str],
        skip_samples: Optional[List[str]] = None,
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
            namespace=namespace or self.namespace,
        )
        self._populate_analysis()
        if ped_files:
            self._populate_pedigree(ped_files)

    def add_stages(self, stages: List[PipelineStage]):
        self.stage_by_name = {
            stage.get_name(): stage for stage in stages
        }

        stage_names = [s.get_name() for s in stages]
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
        
        for i, stage in enumerate(stages):
            if first_stage_num is not None and i < first_stage_num:
                logger.info(f'Skipping stage {stage.get_name()}')
                continue
            logger.info(f'*' * 60)
            logger.info(f'Stage {stage.get_name()}')
            stage.active = True
            stage.add_to_the_pipeline(self)
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


def get_hail_bucket(tmp_bucket, keep_scratch):
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    return hail_bucket


def setup_batch(
    title: str, 
    keep_scratch: bool,
    tmp_bucket: str,
    analysis_project: Project,
    billing_project: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> Batch:
    if not hail_bucket:
        hail_bucket = get_hail_bucket(tmp_bucket, keep_scratch)
    billing_project = (
        billing_project or
        os.getenv('HAIL_BILLING_PROJECT') or
        analysis_project.name
    )
    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
        token=os.getenv('HAIL_TOKEN'),
    )
    b = Batch(name=title, backend=backend)
    return b
