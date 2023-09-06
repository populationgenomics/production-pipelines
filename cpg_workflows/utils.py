"""
Utility functions and constants.
"""

import logging
import hail as hl
import re
import string
import sys
import time
import traceback
import unicodedata
from functools import lru_cache
from os.path import join
from random import choices
from typing import cast, Union

from hailtop.batch import ResourceFile

from cpg_utils import Path, to_path
from cpg_utils.config import get_config


def read_hail(path):
    """
    read a hail object using the appropriate method
    Args:
        path (str): path to the input object
    Returns:
        hail object (hl.MatrixTable or hl.Table)
    """
    if path.strip('/').endswith('.ht'):
        t = hl.read_table(str(path))
    else:
        assert path.strip('/').endswith('.mt')
        t = hl.read_matrix_table(str(path))
    logging.info(f'Read data from {path}')
    return t


def checkpoint_hail(t, file_name: str, checkpoint_prefix: str | None = None):
    if checkpoint_prefix:
        path = join(checkpoint_prefix, file_name)
        if can_reuse(path):
            t = read_hail(str(path))
        else:
            t.write(str(path), overwrite=True)
            logging.info(f'Wrote checkpoint {path}')
            t = read_hail(str(path))
    return t


def missing_from_pre_collected(test: set[Path], known: set[Path]) -> Path | None:
    """
    Check if a path exists in a set of known paths.

    This is useful when checking if a path exists in a set of paths that were
    already collected. This method has been included to permit simple mocking

    Args:
        test (set[Path]): all the files we want to check
        known (set[Path]): all the files we know about

    Returns:
        Path | None: the first path that is missing from the known set, or None
            Path is arbitrary, as the set is unordered
            None indicates No missing files
    """
    return next((p for p in test if p not in known), None)


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


def rich_sequencing_group_id_seds(
    rich_id_map: dict[str, str],
    file_names: list[str | ResourceFile],
) -> str:
    """
    Helper function to add seds into a command that would extend sequencing group IDs
    in each file in `file_names` with an external ID, only if external ID is
    different from the original.

    @param rich_id_map: map used to replace sequencing groups, e.g. {'CPG1': 'CPG1|EXTID'}
    @param file_names: file names and Hail Batch Resource files where to replace IDs
    @return: bash command that does replacement
    """
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f'sed -iBAK \'s/{sgid}/{rich_sgid}/g\' {fname}'
            cmd += '\n'
    return cmd


ExpectedResultT = Union[Path, dict[str, Path], dict[str, str], str, None]
