"""
Extending the Hail's `Batch` class.
"""

import logging
import os
from typing import Optional, Dict

import hailtop.batch as hb
from hailtop.batch.job import Job

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Batch(hb.Batch):
    """
    Thin subclass of the Hail `Batch` class. The aim is to be able to register
    the create jobs to be able to print statistics before submitting.
    """
    def __init__(self, name, backend, *args, **kwargs):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry:
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
    title: str, 
    keep_scratch: bool = False,
    tmp_bucket: Optional[str] = None,
    billing_project: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> Batch:
    """
    Wrapper around the initialization of a Hail Batch object.
    Handles setting the temporary bucket and the billing project.

    :param title: descriptive name of the Batch (will be displayed in the GUI)
    :param billing_project: Hail billing project name
    :param tmp_bucket: path to temporary bucket. Will be used if neither 
        the hail_bucket parameter nor HAIL_BUCKET env var are set.
    :param keep_scratch: whether scratch will be kept after the batch is finished
    :param hail_bucket: bucket for Hail Batch intermediate files.
    """
    if not hail_bucket:
        hail_bucket = get_hail_bucket(tmp_bucket, keep_scratch)

    billing_project = billing_project or os.getenv('HAIL_BILLING_PROJECT')
    if not billing_project:
        raise ValueError(
            'Either the billing_project parameter, or the HAIL_BILLING_PROJECT'
            'environment variable must be set'
        )
    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=hail_bucket,
        token=os.environ.get('HAIL_TOKEN'),
    )
    return Batch(name=title, backend=backend)


def get_hail_bucket(
    tmp_bucket: Optional[str] = None, 
    keep_scratch: bool = False,
) -> str:
    """
    Get bucket where Hail Batch will keep scratch files
    """
    hail_bucket = os.environ.get('HAIL_BUCKET')
    
    if not hail_bucket and not tmp_bucket:
        raise ValueError(
            'Either the tmp_bucket parameter, or the HAIL_BUCKET '
            'environment variable must be set.'
        )

    if keep_scratch and not tmp_bucket:
        raise ValueError(
            'When keep_scratch=True, the tmp_bucket parameter must be set. '
            'Scratch files can be large, so we want to use the tmp bucket '
            'to store them, which is expected to set up to get cleaned '
            'automatically on schedule.'
        )
        
    if keep_scratch or not hail_bucket:
        assert tmp_bucket
        hail_bucket = os.path.join(tmp_bucket, 'hail')

    return hail_bucket


def job_name(name, sample: str = None, dataset: str = None) -> str:
    """
    Extend the descriptive job name to reflect the dataset and the sample names
    """
    if sample and dataset:
        name = f'{dataset}/{sample}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    return name
