import os
import shutil
import tempfile
from dataclasses import dataclass
from enum import Enum
import logging
from typing import List, Dict, Optional, Any, Tuple, Callable, Union
from abc import ABC, abstractmethod

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines import utils
from cpg_production_pipelines.smdb import SMDB, Analysis, Sequence
from cpg_production_pipelines.hailbatch import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Namespace(Enum):
    TMP = 1
    TEST = 2
    MAIN = 3


class AnalysisType(Enum):
    """Types of analysis"""

    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'


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
            sample = f'{project}/{sample}'
        if sample:
            name = f'{sample}: {name}'

        if label and (sample or samples):
            if label not in self.labelled_jobs:
                self.labelled_jobs[label] = {'job_n': 0, 'samples': set()}
            self.labelled_jobs[label]['job_n'] += 1
            self.labelled_jobs[label]['samples'] |= (samples or {sample})
        else:
            self.other_job_num += 1
        self.total_job_num += 1
        j = super().new_job(name, attributes=attributes, **kwargs)
        return j
    
    
# @dataclass
# class StageOutput:
#     expected_path: str
#     found_path: str


class PipelineStage(ABC):
    def __init__(
        self, 
        pipe: 'Pipeline', 
        name: str,
        analysis_type: AnalysisType = AnalysisType.CUSTOM,
    ):
        self.pipe = pipe
        self.name = name
        self.analysis_type = analysis_type,
        self.jobs: List[Job] = list()
        self.inputs: Dict[str, str] = dict()
        self.outputs: Dict[str, str] = dict()
        # self._jobs: Dict[Tuple[str, str, str], Job] = dict()
        # self._expected_outputs: Dict[Tuple[str, str, str], str] = dict()
        # self._found_outputs: Dict[Tuple[str, str, str], str] = dict()

    # def get_output(self, name='', sid='', pid='') -> Optional[str]:
    #     return self._outputs.get((name, sid, pid))
    # 
    # def _set_output(self, ouptut_path: str, name='', sid='', pid=''):
    #     self._outputs[(name, sid, pid)] = ouptut_path
    # 
    # def get_job(self, name='', sid='', pid='') -> Optional[Job]:
    #     return self._jobs.get((name, sid, pid))
    # 
    # def _set_job(self, job: Job, name='', sid='', pid=''):
    #     self._jobs[(name, sid, pid)] = job
    # 
    # def get_all_jobs(self) -> List[Job]:
    #     return list(self._jobs.values())

    @abstractmethod
    def add_jobs(self, **kwargs):
        pass

    @abstractmethod
    def get_expected_output(self, *args):
        pass
    
    def find_stage_output(
        self, obj: Union['Sample', 'Project', 'Pipeline'],
    ) -> Optional[str]:
        found_path = obj.output_by_stage.get(self.name)
        if not found_path:
            analysis = obj.analysis_by_type.get(self.analysis_type)
            if not analysis:
                return None

            if self.pipe.validate_smdb_analyses:
                if isinstance(obj, 'Sample'):
                    sample_ids = [obj.id]
                elif isinstance(obj, 'Project'):
                    sample_ids = [s.id for s in obj.samples]
                else:
                    sample_ids = obj.get_all_sample_ids()
                
                expected_path = self.get_expected_output(obj)
                found_path = self.pipe.db.process_existing_analysis(
                    sample_ids=sample_ids,
                    completed_analysis=analysis,
                    analysis_type=self.analysis_type,
                    expected_output_fpath=expected_path,
                )
            else:
                found_path = analysis.output
            obj.output_by_stage[self.name] = found_path
        return found_path
    
    
class SampleStage(PipelineStage, ABC):
    """
    Sample-level stage. run() takes sample and project name as input
    """
    @abstractmethod
    def add_jobs(self, sample: 'Sample'):
        pass


class ProjectStage(PipelineStage, ABC):
    """
    Project-level stage. run() takes project name as input
    """
    @abstractmethod
    def add_jobs(self, project: 'Project'):
        pass


class CohortStage(PipelineStage, ABC):    
    """
    Entire cohort level stage
    """
    @abstractmethod
    def add_jobs(self):
        pass


@dataclass
class PedigreeInfo:
    fam_id: str
    dad: Optional['Sample']
    mom: Optional['Sample']
    sex: str
    phenotype: str


@dataclass
class Sample:
    """
    Represents a Sample.
    """
    id: str
    external_id: str
    project: 'Project'
    alignment_input: Optional[AlignmentInput] = None
    good: bool = True
    seq_info: Optional[Sequence] = None
    # From SMDB Analysis entries
    analysis_by_type: Dict[AnalysisType, Analysis] = None
    # Other outputs, not supported by DB.
    # Project-level outputs are specified in Project class,
    # Cohort-level outputs are specified in Pipeline class
    output_by_stage: Dict[str, str] = None
    jobs_by_stage: Dict[str, List[Job]] = None
    pedigree: Optional[PedigreeInfo] = None 


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


@dataclass
class Project:
    """
    Represents a CPG dataset/project.
    """
    name: str
    samples: List[Sample]
    pipeline: 'Pipeline'
    is_test: bool = False
    # From SMDB Analysis entries
    analysis_by_type: Dict[AnalysisType, Analysis] = None
    # Outputs not supported by DB
    # Sample-level outputs are specified and Sample class,
    # Cohort-level outputs are specified in Pipeline class
    output_by_stage: Dict[str, str] = None
    jobs_by_stage: Dict[str, Job] = None

    def get_bucket(self):
        return f'gs://cpg-{self.name}-{self.pipeline.output_suf}'

    def get_tmp_bucket(self):
        return f'gs://cpg-{self.name}-{self.pipeline.output_suf}-tmp'


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
        smdb_update_analyses: bool = False,
        smdb_check_existence: bool = False,
        validate_smdb_analyses: bool = False,
        hail_billing_project: Optional[str] = None,
        overwrite: bool = True,
        first_stage: Optional[str] = None,
        last_stage: Optional[str] = None,
        config: Any = None,
    ):
        self.analysis_project = analysis_project
        self.name = name
        self.output_version = output_version
        self.namespace = namespace
        self.overwrite = overwrite
        self.first_stage = first_stage
        self.last_stage = last_stage
        self.db = SMDB(
            self.analysis_project,
            do_update_analyses=smdb_update_analyses,
            do_check_existence=smdb_check_existence,
        )
        self.validate_smdb_analyses = validate_smdb_analyses
        self.projects: List[Project] = []
        self.local_tmp_dir = tempfile.mkdtemp()
        self.keep_scratch = keep_scratch
        self.config = config

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
            f'gs://cpg-{self.analysis_project}-{{suffix}}/'
            f'{self.name}/'
            f'{self.output_version}'
        )
        self.tmp_bucket = path_ptrn.format(suffix=tmp_suf)
        self.analysis_bucket = path_ptrn.format(suffix=analysis_suf)
        self.web_bucket = path_ptrn.format(suffix=web_suf)

        self.b = setup_batch(
            title, 
            self.keep_scratch, 
            self.tmp_bucket, 
            self.analysis_project, 
            hail_billing_project
        )

        # Found global SMDB analysis entries by type
        self.analysis_by_type: Dict[AnalysisType, Analysis] = None
        # Cohort-level outputs, not supported by SMDB.
        # Sample- and project-level outputs are specified in Sample and Project classes
        self.output_by_stage: Dict[str, str] = None
        self.jobs_by_stage: Dict[str, Job] = None

        self.stages = List[PipelineStage]

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
        skip_samples: Optional[List[str]] = None,
        namespace: Optional[str] = None,
    ):
        namespace = namespace or self.output_suf
        samples_by_project = self.db.get_samples_by_project(
            project_names=input_projects,
            namespace=namespace,
            skip_samples=skip_samples,
        )
        for proj_name, sample_datas in samples_by_project.items():
            project = Project(
                name=proj_name,
                samples=[],
                pipeline=self,
                is_test=namespace != 'main',
                output_by_stage=dict(),
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
        ),

        for project in self.projects:
            sample_ids = [s.id for s in project.samples]
            
            seq_info_by_sid = self.db.find_seq_info_by_sid(sample_ids)

            cram_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='cram',
            )
            gvcf_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='gvcf',
            )

            for s in project.samples:
                s.seq_info = seq_info_by_sid.get(s.id)
                s.analysis_by_type[AnalysisType.CRAM] = cram_per_sid.get(s.id),
                s.analysis_by_type[AnalysisType.GVCF] = gvcf_per_sid.get(s.id),

        self.analysis_by_type[AnalysisType.JOINT_CALLING] = jc_analysis
            
    def _populate_pedigree(self, ped_files: List[str]):
        sample_by_extid = dict()
        for s in self.get_all_samples():
            sample_by_extid[s.external_id] = s
        
        for ped_file in ped_files:
            with open(ped_file) as f:
                for line in f:
                    if 'Family.ID' in line:
                        continue
                    fields = line.strip().split()[:6]
                    fam_id, sam_id, pat_id, mat_id, sex, phenotype = fields
                    if sam_id in sample_by_extid:
                        s = sample_by_extid[sam_id]
                        s.pedigree = PedigreeInfo(
                            fam_id=fam_id,
                            dad=sample_by_extid.get(pat_id),
                            mom=sample_by_extid.get(mat_id),
                            sex=sex,
                            phenotype=phenotype,
                        )
        for project in self.projects:
            samples_with_ped = [s for s in project.samples if s.pedigree]
            logger.info(
                f'{project}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(project.samples)}'
            )

    def populate_samples(
        self,
        input_projects: List[str],
        skip_samples: Optional[List[str]] = None,
        namespace: Optional[str] = None,
        ped_files: Optional[List[str]] = None,
    ) -> None:
        """
        Finds input samples, analyses and sequences from the DB, 
        populates self.projects, adds pedigree information
        """
        self._populate_projects(
            input_projects=input_projects,
            skip_samples=skip_samples,
            namespace=namespace
        )
        self._populate_analysis()
        if ped_files:
            self._populate_pedigree(ped_files)

    def add_sample_stage(self, stage: SampleStage):
        for project in self.projects:
            logger.info(f'{stage.name}: adding jobs for project {project.name}')
            for sample in project.samples:
                logger.info(f'{stage.name}: adding jobs for {project.name}/{sample.id}')
                stage.add_jobs(sample)

    def add_project_stage(self, stage: ProjectStage):
        for project in self.projects:
            logger.info(f'{stage.name}: adding jobs for project {project.name}')
            stage.add_jobs(project)

    def add_cohort_stage(self, stage: CohortStage):
        stage.add_jobs()

    def add_stages(self, stage_by_name: Dict[str, PipelineStage]):
        for name, stage in stage_by_name.items():
            if name == self.config.last_stage:
                logger.info(f'Last stage is {name}, stopping here')
            elif self.config.first_stage and name != self.config.first_stage:
                logger.info(f'Skipping stage {name}')
                self.config.first_stage = None
            else:
                logger.info(f'*' * 60)
                logger.info(f'Stage {name}')
                if isinstance(stage, SampleStage):
                    self.add_sample_stage(stage)
                elif isinstance(stage, ProjectStage):
                    self.add_project_stage(stage)
                else:
                    self.add_cohort_stage(stage)
                logger.info(f'')
    
    def can_reuse(self, fpath: Optional[str]) -> bool:
        """
        Checks if the fpath exists, but always returns False if self.overwrite=True
        """
        return utils.can_reuse(fpath, overwrite=self.overwrite)


def setup_batch(
    title, 
    keep_scratch,
    tmp_bucket,
    analysis_project,
    billing_project: Optional[str] = None,
) -> Batch:
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    billing_project = (
        billing_project or
        os.getenv('HAIL_BILLING_PROJECT') or
        analysis_project
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
