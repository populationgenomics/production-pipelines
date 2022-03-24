"""
Extending the Hail's `Batch` class.
"""

import logging
import os

import hailtop.batch as hb
from cloudpathlib import CloudPath
from hailtop.batch.job import Job

from .. import Path, to_path

logger = logging.getLogger(__file__)


class Batch(hb.Batch):
    """
    Thin subclass of the Hail `Batch` class. The aim is to be able to register
    created jobs, in order to print statistics before submitting the Batch.
    """
    def __init__(self, name, backend, *args, **kwargs):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry:
        self.labelled_jobs = dict()
        self.other_job_num = 0
        self.total_job_num = 0

    def new_job(
        self,
        name: str | None = None,
        attributes: dict[str, str] | None = None,
        **kwargs,
    ) -> Job:
        """
        Adds job to the Batch, and also registers it in `self.job_stats` for
        statistics.
        """
        if not name:
            logger.critical('Error: job name must be defined')
        
        attributes = attributes or dict()
        dataset = attributes.get('dataset')
        sample = attributes.get('sample')
        samples = attributes.get('samples')
        label = attributes.get('label', name)

        name = job_name(name, sample, dataset)

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


def setup_batch(
    description: str, 
    billing_project: str | None = None,
    hail_bucket: Path | None = None,
) -> Batch:
    """
    Wrapper around the initialization of a Hail Batch object.
    Handles setting the temporary bucket and the billing project.

    @param description: descriptive name of the Batch (will be displayed in the GUI)
    @param billing_project: Hail billing project name
    @param hail_bucket: bucket for Hail Batch intermediate files.
    """
    billing_project = get_billing_project(billing_project)
    hail_bucket = get_hail_bucket(hail_bucket)

    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=str(hail_bucket),
        token=os.environ.get('HAIL_TOKEN'),
    )
    return Batch(name=description, backend=backend)


def get_hail_bucket(hail_bucket: str | Path | None = None) -> Path:
    """
    Get Hail bucket.
    """
    hail_bucket = hail_bucket or os.getenv('HAIL_BUCKET')
    if not hail_bucket:
        raise ValueError(
            'Either the hail_bucket parameter, or the HAIL_BUCKET'
            'environment variable must be set'
        )
    return to_path(hail_bucket)


def get_billing_project(billing_project: str | None = None) -> str:
    """
    Get Hail billing project.
    """
    
    billing_project = billing_project or os.getenv('HAIL_BILLING_PROJECT')
    if not billing_project:
        raise ValueError(
            'Either the billing_project parameter, or the HAIL_BILLING_PROJECT'
            'environment variable must be set'
        )
    return billing_project


def job_name(name, sample: str = None, dataset: str = None) -> str:
    """
    Extend the descriptive job name to reflect the dataset and the sample names
    """
    if sample and dataset:
        name = f'{dataset}/{sample}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    return name


def hail_query_env(j: Job, hail_billing_project: str, hail_bucket: Path | None = None):
    """
    Setup environment to run Hail Query Service backend script.
    """
    j.env('HAIL_BILLING_PROJECT', hail_billing_project)
    hail_bucket = get_hail_bucket(hail_bucket)
    # Has to be without prefix:
    if isinstance(hail_bucket, CloudPath):
        bucket = str(hail_bucket)[len(hail_bucket.cloud_prefix):]
    else:
        bucket = str(hail_bucket)
    j.env('HAIL_BUCKET', bucket)
    j.env('HAIL_SHA', os.environ['HAIL_SHA'])
    j.env('HAIL_JAR_URL', os.environ['HAIL_JAR_URL'])
