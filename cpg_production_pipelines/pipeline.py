import os
import shutil
import tempfile
from dataclasses import dataclass
from enum import Enum
import logging
from typing import List, Dict, Optional, Set

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines.smdb import SMDB
from cpg_production_pipelines.utils import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Namespace(Enum):
    TMP = 1
    TEST = 2
    MAIN = 3


@dataclass
class Sample:
    """
    Represents a Sample.
    """
    id: str
    external_id: str
    project: str
    alignment_input: Optional[AlignmentInput] = None
    good: bool = True
   
    
@dataclass
class Project:
    """
    Represents a CPG dataset/project.
    """
    name: str
    samples: List[Sample]
    is_test: bool = False


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
        smdb_check_existence: bool = False
    ):
        self.analysis_project = analysis_project
        self.name = name
        self.output_version = output_version
        self.namespace = namespace
        self.db = SMDB(
            self.analysis_project,
            do_update_analyses=smdb_update_analyses,
            do_check_existence=smdb_check_existence,
        )

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

        self.keep_scratch = keep_scratch
        self.b = setup_batch(self.tmp_bucket, title, keep_scratch)
        self.local_tmp_dir = tempfile.mkdtemp()

        self.projects: List[Project] = []
        
        # Job stats registry
        self.labelled_jobs = dict()
        self.other_job_num = 0
        self.total_job_num = 0

    def get_all_samples(self) -> List[Sample]:
        all_samples = []
        for proj in self.projects:
            all_samples.extend(proj.samples)
        return all_samples

    def add_job(
        self, 
        name: str, 
        label: Optional[str], 
        samples: Optional[Set[Sample]],
    ) -> Job:
        """
        Adds job to the Batch, and also registers it in `self.job_stats` for
        statistics.
        """
        if label and samples is not None:
            if label not in self.labelled_jobs:
                self.labelled_jobs[label] = {'job_n': 0, 'samples': set()}
            self.labelled_jobs[label]['job_n'] += 1
            self.labelled_jobs[label]['samples'] |= samples
        else:
            self.other_job_num += 1
        self.total_job_num += 1
        return self.b.new_job(name)

    def run(self, dry_run: bool = False) -> None:
        if self.b:
            logger.info(f'Will submit {self.total_job_num} jobs:')
            for label, stat in self.labelled_jobs.items():
                logger.info(f'  {label}: {stat["job_n"]} for {len(stat["samples"])} samples')
            logger.info(f'  Other jobs: {self.other_job_num}')

            self.b.run(
                dry_run=dry_run,
                delete_scratch_on_exit=not self.keep_scratch,
                wait=False,
            )
        shutil.rmtree(self.local_tmp_dir)

    def find_samples(
        self,
        input_projects: List[str],
        skip_samples: Optional[List[str]] = None,
        namespace: Optional[str] = None,
    ) -> None:
        """
        Populates self.projects
        """
        namespace = namespace or self.output_suf
        samples_by_project = self.db.get_samples_by_project(
            projects=input_projects,
            namespace=namespace,
            skip_samples=skip_samples,
        )
        for proj_name, samples in samples_by_project.items():
            p = Project(
                proj_name, 
                is_test=namespace != 'main',
                samples=[Sample(
                    id=s['id'], 
                    external_id=s['external_id'],
                    project=proj_name,
                ) for s in samples]
            )
            self.projects.append(p)


def setup_batch(tmp_bucket, title, keep_scratch) -> hb.Batch:
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
        token=os.getenv('HAIL_TOKEN'),
    )
    b = hb.Batch(title, backend=backend)
    return b
