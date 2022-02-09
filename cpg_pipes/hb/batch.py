"""
Extending the Hail's `Batch` class.
"""

import logging
import os
from typing import Optional

import hailtop.batch as hb

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def setup_batch(
    title: str, 
    keep_scratch: bool = False,
    tmp_bucket: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> hb.Batch:
    """
    Wrapper around the initialization of a Hail Batch object.
    Handles setting the temporary bucket and the billing project.

    :param title: descriptive name of the Batch (will be displayed in the GUI)
    :param tmp_bucket: path to temporary bucket. Will be used if neither 
        the hail_bucket parameter nor HAIL_BUCKET env var are set.
    :param keep_scratch: whether scratch will be kept after the batch is finished
    :param hail_bucket: bucket for Hail Batch intermediate files.
    """
    if not hail_bucket:
        hail_bucket = get_hail_bucket(tmp_bucket, keep_scratch)

    billing_project = os.getenv('HAIL_BILLING_PROJECT')
    if not billing_project:
        raise ValueError(
            'HAIL_BILLING_PROJECT environment variable must be set'
        )
    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
        token=os.environ.get('HAIL_TOKEN'),
    )
    return hb.Batch(name=title, backend=backend)


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


def _job_name(name, sample: str = None, project: str = None) -> str:
    """
    Extend the descriptive job name to reflect the project and the sample names
    """
    if sample and project:
        name = f'{project}/{sample}: {name}'
    elif project:
        name = f'{project}: {name}'
    return name
