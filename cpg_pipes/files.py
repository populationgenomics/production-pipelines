"""
Checking existence of intermediate files that can be on buckets.
"""

import logging
import os
from typing import Optional, Union, Iterable
from google.cloud import storage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def exists(path: str, silent: bool = False) -> bool:
    """
    Check if the object exists, where the object can be:
    - local file
    - local directory
    - cloud object
    - cloud path to a *.mt or *.ht Hail Query table,
      in which case it will check for the existence of a
      *.mt/_SUCCESS or *.ht/_SUCCESS file
    For cloud objects, the function would log each object access, which is useful
    to detect unnecessary accesses which can be a bottleneck when working with a lot
    of files. The `silent` flag can be set to True to disable logging.
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        p = path.replace('gs://', '').split('/', maxsplit=1)[1]
        p = p.rstrip('/')  # ".mt/" -> ".mt"
        if any(p.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            p = os.path.join(p, '_SUCCESS')
        gs = storage.Client()
        exists = gs.get_bucket(bucket).get_blob(p) is not None
        if silent:
            if exists:
                logger.info(f'Checking object existence, exists: {path}')
            else:
                logger.info(f'Checking object existence, doesn\'t exist: {path}')
        return exists
    return os.path.exists(path)


def can_reuse(
    path: Optional[Union[Iterable[str], str]],
    overwrite: bool,
    silent=False,
) -> bool:
    """
    Checks if `path` is good to reuse in a pipelines: the object exists
    and `overwrite` is False.

    If `path` is a collection, it requires all files in it to exist.
    """
    if overwrite:
        return False

    if not path:
        return False

    if not isinstance(path, str):
        return all(can_reuse(fp, overwrite) for fp in path)

    if not exists(path, silent=silent):
        return False

    if not silent:
        logger.info(f'Reusing existing {path}. Use --overwrite to overwrite')
    return True
