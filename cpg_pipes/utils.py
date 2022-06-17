"""
Utility functions and constants.
"""

import logging
import string
import sys
import time
import traceback
from functools import lru_cache
from random import choices
from typing import cast

from . import Path, to_path

logger = logging.getLogger(__file__)

# Packages to install on a dataproc cluster, to use with the dataproc wrapper.
DATAPROC_PACKAGES = [
    'cpg_utils==4.3.6.2',
    'cpg_pipes',
    'cpg_gnomad',  # github.com/populationgenomics/gnomad_methods
    'seqr_loader==1.2.5',  # hail-elasticsearch-pipelines
    'elasticsearch==8.1.1',
    'cpg_utils',
    'click',
    'google',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]


@lru_cache
def exists(path: Path | str, verbose: bool = True) -> bool:
    """
    Caching version of the existence check. 
    The python code runtime happens entirely during the pipeline submittion, 
    without waiting for it to finish, so there is no expectation that object 
    existence status would change during the runtime. This, this function uses
    `@lru_cache` to make sure that object existence is checked only once.
    """
    return exists_not_cached(path, verbose)


def exists_not_cached(path: Path | str, verbose: bool = True) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * cloud object
        * cloud URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    @param path: path to the file/directory/object/mt/ht
    @param verbose: print on each check
    @return: True if the object exists
    """
    path = cast(Path, to_path(path))

    # rstrip to ".mt/" -> ".mt"
    if any(str(path).rstrip('/').endswith(f'.{suf}') for suf in ['mt', 'ht']):
        path = path / '_SUCCESS'

    if verbose:
        # noinspection PyBroadException
        try:
            res = path.exists()
        except BaseException:
            traceback.print_exc()
            logger.error(f'Failed checking {path}')
            sys.exit(1)
        logger.debug(f'Checked {path} [' + ('exists' if res else 'missing') + ']')
        return res
    return path.exists()


def can_reuse(
    path: list[Path] | Path | str | None,
    overwrite: bool,
) -> bool:
    """
    Checks if `fpath` is good to reuse in the analysis: it exists
    and `overwrite` is False.

    If `fpath` is a collection, it requires all files in it to exist.
    """
    if overwrite:
        return False

    if not path:
        return False

    if isinstance(path, list):
        return all(can_reuse(fp, overwrite) for fp in path)

    if not exists(path):
        return False

    logger.debug(f'Reusing existing {path}. Use --overwrite to overwrite')
    return True


def timestamp(rand_suffix_len: int = 5) -> str:
    """
    Generate a timestamp string. If `rand_suffix_len` is set, adds a short random 
    string of this length for uniqueness.
    """
    result = time.strftime('%Y_%m%d_%H%M')
    if rand_suffix_len:
        rand_bit = ''.join(choices(
            string.ascii_uppercase + string.digits, k=rand_suffix_len)
        )
        result += f'_{rand_bit}' 
    return result
