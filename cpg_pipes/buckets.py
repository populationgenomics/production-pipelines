"""
Utility functions to interact with objects on buckets.
"""

import logging
import os
import subprocess
from pathlib import Path
from typing import Iterable, cast

from google.cloud import storage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def file_exists(path: str | Path) -> bool:
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
    if isinstance(path, Path):
        path = cast(str, path.absolute())

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
    fpath: Iterable[str|Path]|str|Path|None,
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

    if isinstance(fpath, Iterable):
        return all(can_reuse(fp, overwrite) for fp in fpath)

    if isinstance(fpath, Path) and not fpath.exists():
        return False
    
    if isinstance(fpath, str) and not file_exists(fpath):
        return False

    if not silent:
        logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
    return True


def gsutil_cp(
    src_path: str,
    dst_path: str,
    disable_check_hashes: bool = False,
    recursive: bool = False,
    quiet: bool = False,
) -> str:
    """
    Wrapper around `gsutil cp`

    :param src_path: path to a file to copy from
    :param dst_path: path to copy to
    :param disable_check_hashes:
        Uses the gsutil option `-o GSUtil:check_hashes=never` which is required to
        get around the gsutil integrity checking error, as conda gsutil doesn't use
        CRC32c:
        > Downloading this composite object requires integrity checking with CRC32c,
          but your crcmod installation isn't using the module's C extension, so the
          hash computation will likely throttle download performance.

          To download regardless of crcmod performance or to skip slow integrity
          checks, see the "check_hashes" option in your boto config file.
    :param recursive: to copy a directory
    :param quiet: disable logging of commands and copied files
    :returns: dst_path
    """
    cmd = (
        'gsutil '
        + ('-q ' if quiet else '')
        + ('-o GSUtil:check_hashes=never ' if disable_check_hashes else '')
        + 'cp '
        + ('-r ' if recursive else '')
        + f'{src_path} {dst_path}'
    )
    if not quiet:
        logger.info(cmd)
    subprocess.run(cmd, check=False, shell=True)
    return dst_path
