import os
import shutil
import tempfile
from dataclasses import dataclass
from enum import Enum
import logging
from typing import List, Dict

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines.smdb import SMDB

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Namespace(Enum):
    TMP = 1
    TEST = 2
    MAIN = 3


@dataclass
class Sample:
    id: str
    external_id: str
    project: str


class Pipeline:
    def __init__(
        self,
        analysis_project: str,
        name: str,
        output_version: str,
        namespace: Namespace,
        keep_scratch: bool,
        title: str,
        check_inputs_existence: bool = True,
    ):
        self.analysis_project = analysis_project
        self.name = name
        self.output_version = output_version
        self.namespace = namespace

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

        self.samples_by_project: Dict[str, Sample] = dict()
        self.check_inputs_existence: bool = check_inputs_existence

    def get_all_samples(self) -> List[Sample]:
        all_samples = []
        for proj, samples in self.samples_by_project.items():
            all_samples.extend(samples)
        return all_samples

    def add_job(self, name: str) -> Job:
        return self.b.new_job(name)

    def run(self, dry_run):
        if self.b:
            self.b.run(
                dry_run=dry_run,
                delete_scratch_on_exit=not self.keep_scratch,
                wait=False,
            )
        shutil.rmtree(self.local_tmp_dir)


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
