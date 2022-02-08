"""
Utility functions to interact with objects on buckets.
"""

import logging
import os
from typing import Optional, Union, Iterable

from google.cloud import storage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
        * Google Storage URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    :param path: path to the file/directory/object/mt/ht
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        p = path.replace('gs://', '').split('/', maxsplit=1)[1]
        p = p.rstrip('/')  # ".mt/" -> ".mt"
        if any(p.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            p = os.path.join(p, '_SUCCESS')
        gs = storage.Client()
        exists = gs.get_bucket(bucket).get_blob(p) is not None
        if exists:
            logger.info(f'Checking object existence, exists: {path}')
        else:
            logger.info(f'Checking object existence, doesn\'t exist: {path}')
        return exists
    return os.path.exists(path)


def can_reuse(
    fpath: Optional[Union[Iterable[str], str]],
    overwrite: bool,
    silent=False,
) -> bool:
    """
    Checks if `fpath` is good to reuse in the analysis: it exists
    and `overwrite` is False.

    If `fpath` is a collection, it requires all files in it to exist.
    """
    if overwrite:
        return False

    if not fpath:
        return False

    if not isinstance(fpath, str):
        return all(can_reuse(fp, overwrite) for fp in fpath)

    if not file_exists(fpath):
        return False

    if not silent:
        logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
    return True
