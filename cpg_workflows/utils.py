"""
Utility functions and constants.
"""

import logging
import math
import re
import string
import sys
import time
import traceback
import unicodedata
from functools import lru_cache
from random import choices
from typing import cast

from hailtop.batch import ResourceFile

from cpg_utils import Path, to_path
from cpg_utils.config import get_config


@lru_cache
def exists(path: Path | str, verbose: bool = True) -> bool:
    """
    `exists_not_cached` that caches the result.

    The python code runtime happens entirely during the workflow construction,
    without waiting for it to finish, so there is no expectation that the object
    existence status would change during the runtime. This, this function uses
    `@lru_cache` to make sure that object existence is checked only once.
    """
    return exists_not_cached(path, verbose)


def exists_not_cached(path: Path | str, verbose: bool = True) -> bool:
    """
    Check if the object by path exists, where the object can be:
        * local file,
        * local directory,
        * cloud object,
        * cloud or local *.mt, *.ht, or *.vds Hail data, in which case it will check
          for the existence of a corresponding _SUCCESS object instead.
    @param path: path to the file/directory/object/mt/ht
    @param verbose: print on each check
    @return: True if the object exists
    """
    path = cast(Path, to_path(path))

    if path.suffix in ['.mt', '.ht']:
        path /= '_SUCCESS'
    if path.suffix in ['.vds']:
        path /= 'variant_data/_SUCCESS'

    if verbose:
        # noinspection PyBroadException
        try:
            res = path.exists()
        except BaseException:
            traceback.print_exc()
            logging.error(f'Failed checking {path}')
            sys.exit(1)
        logging.debug(f'Checked {path} [' + ('exists' if res else 'missing') + ']')
        return res
    return path.exists()


def can_reuse(
    path: list[Path] | Path | str | None,
    overwrite: bool = False,
) -> bool:
    """
    Checks if the object at `path` is good to reuse:
    * overwrite has the default value of False,
    * check_intermediates has the default value of True,
    * object exists.

    If `path` is a collection, it requires all paths to exist.
    """
    if overwrite:
        return False

    if not get_config()['workflow'].get('check_intermediates', True):
        return False

    if not path:
        return False

    paths = path if isinstance(path, list) else [path]
    if not all(exists(fp, overwrite) for fp in paths):
        return False

    logging.debug(f'Reusing existing {path}')
    return True


def timestamp(rand_suffix_len: int = 5) -> str:
    """
    Generate a timestamp string. If `rand_suffix_len` is set, adds a short random
    string of this length for uniqueness.
    """
    result = time.strftime('%Y_%m%d_%H%M')
    if rand_suffix_len:
        rand_bit = ''.join(
            choices(string.ascii_uppercase + string.digits, k=rand_suffix_len)
        )
        result += f'_{rand_bit}'
    return result


def slugify(line: str):
    """
    Slugify a string.

    Example:
    >>> slugify(u'Héllø W.1')
    'hello-w-1'
    """

    line = unicodedata.normalize('NFKD', line).encode('ascii', 'ignore').decode()
    line = line.strip().lower()
    line = re.sub(
        r'[\s.]+',
        '-',
        line,
    )
    return line


def rich_sample_id_seds(
    rich_id_map: dict[str, str],
    file_names: list[str | ResourceFile],
) -> str:
    """
    Helper function to add seds into a command that would extend samples IDs
    in each file in `file_names` with an external ID, only if external ID is
    different from the original.

    @param rich_id_map: map used to replace samples, e.g. {'CPG1': 'CPG1|EXTID'}
    @param file_names: file names and Hail Batch Resource files where to replace IDs
    @return: bash command that does replacement
    """
    cmd = ''
    for sid, rich_sid in rich_id_map.items():
        for fname in file_names:
            cmd += f'sed -iBAK \'s/{sid}/{rich_sid}/g\' {fname}'
            cmd += '\n'
    return cmd
