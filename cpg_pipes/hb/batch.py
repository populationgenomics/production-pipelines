"""
Extending the Hail's `Batch` class.
"""

import logging
import os

from cloudpathlib import CloudPath
import hailtop.batch as hb
from hailtop.batch.job import Job, PythonJob

from .. import Path, to_path
from ..providers import Cloud

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
    
    def _process_attributes(
        self,
        name: str | None = None,
        attributes: dict[str, str] | None = None,
    ):
        """
        Use job attributes to make the job name more descriptibe, and add
        labels for Batch pre-submission stats.
        """
        if not name:
            raise ValueError('Error: job name must be defined')

        attributes = attributes or dict()
        dataset = attributes.get('dataset')
        sample = attributes.get('sample')
        samples = attributes.get('samples')
        intervals = attributes.get('intervals')
        label = attributes.get('label', name)

        name = make_job_name(name, sample, dataset, intervals)

        if label and (sample or samples):
            if label not in self.labelled_jobs:
                self.labelled_jobs[label] = {'job_n': 0, 'samples': set()}
            self.labelled_jobs[label]['job_n'] += 1
            self.labelled_jobs[label]['samples'] |= (samples or {sample})
        else:
            self.other_job_num += 1
        self.total_job_num += 1
        return name

    def new_python_job(
        self,
        name: str | None = None,
        attributes: dict[str, str] | None = None,
    ) -> PythonJob:
        """
        Wrapper around `new_python_job()` that processes job attributes.
        """
        name = self._process_attributes(name, attributes)
        return super().new_python_job(name, attributes=attributes)

    def new_job(
        self,
        name: str | None = None,
        attributes: dict[str, str] | None = None,
        **kwargs,
    ) -> Job:
        """
        Wrapper around `new_job()` that processes job attributes.
        """
        name = self._process_attributes(name, attributes)
        return super().new_job(name, attributes=attributes)


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
        remote_tmpdir=hail_bucket,
        token=os.environ.get('HAIL_TOKEN'),
    )
    return Batch(name=description, backend=backend)


def get_hail_bucket(hail_bucket: str | Path | None = None, cloud=Cloud.GS) -> str:
    """
    Get Hail bucket.
    """
    if not hail_bucket:
        hail_bucket = os.getenv('HAIL_BUCKET')
        if not hail_bucket:
            raise ValueError(
                'Either the hail_bucket parameter, or the HAIL_BUCKET '
                'environment variable must be set'
            )
    hail_bucket = str(hail_bucket)
    prefs = {Cloud.GS: 'gs', Cloud.AZ: 'hail-az'}
    if not any(hail_bucket.startswith(f'{pref}://') for pref in prefs.values()):
        hail_bucket = f'{prefs[cloud]}://{hail_bucket}'

    return hail_bucket


def get_billing_project(billing_project: str | None = None) -> str:
    """
    Get Hail billing project.
    """
    
    billing_project = billing_project or os.getenv('HAIL_BILLING_PROJECT')
    if not billing_project:
        raise ValueError(
            'Either the billing_project parameter, or the HAIL_BILLING_PROJECT '
            'environment variable must be set'
        )
    return billing_project


def make_job_name(
    name: str, 
    sample: str | None = None, 
    dataset: str | None = None,
    intervals: str | None = None,
) -> str:
    """
    Extend the descriptive job name to reflect job attributes.
    """
    if sample and dataset:
        name = f'{dataset}/{sample}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    if intervals:
        name += f', {intervals}'
    return name


def hail_query_env(
    j: Job, 
    hail_billing_project: str, 
    hail_bucket: str | Path | None = None
):
    """
    Setup environment to run Hail Query Service backend script.
    """
    j.env('HAIL_BILLING_PROJECT', hail_billing_project)
    j.env('HAIL_BUCKET', get_hail_bucket(hail_bucket))
