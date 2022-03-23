"""
Utility functions and constants.
"""

import hashlib
import logging
import os
import sys
import time
import click
from typing import Any, Callable, cast

from . import Path, to_path, get_package_path

logger = logging.getLogger(__file__)

# Default reference genome build.
DEFAULT_REF = 'GRCh38'

# Packages to install on a dataproc cluster, to use with the dataproc wrapper.
DATAPROC_PACKAGES = [
    # 'cpg_pipes==0.3.0',
    'cpg_gnomad',   # github.com/populationgenomics/gnomad_methods
    'seqr_loader==1.2.0',  # hail-elasticsearch-pipelines
    'elasticsearch==8.1.1',
    'cpg_utils',
    'click',
    'google',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]

# Location of python scripts to be called directly from command line.
SCRIPTS_DIR = to_path('scripts')

# Location of Hail Query scripts, to use with the dataproc wrapper.
QUERY_SCRIPTS_DIR = to_path('query_scripts')

# This python package name.
PACKAGE_DIR = get_package_path()


def get_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    accompanying_metadata_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    @param ext: check that the path has the expected extension
    @param must_exist: check that the input file/object/directory exists
    @param accompanying_metadata_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. genomes.mt and genomes.metadata.ht)
    @return: a callback suitable for Click parameter initialization
    """
    def _callback(_: click.Context, param: click.Option, value: Any):
        if value is None:
            return value
        if ext:
            assert isinstance(value, str), value
            value = value.rstrip('/')
            if not value.endswith(f'.{ext}'):
                raise click.BadParameter(
                    f'The argument {param.name} is expected to have '
                    f'an extension .{ext}, got: {value}'
                )
        if must_exist:
            if not exists(value):
                raise click.BadParameter(f"{value} doesn't exist or incomplete")
            if accompanying_metadata_suffix:
                accompanying_metadata_fpath = (
                    os.path.splitext(value)[0] + accompanying_metadata_suffix
                )
                if not exists(accompanying_metadata_fpath):
                    raise click.BadParameter(
                        f"An accompanying file {accompanying_metadata_fpath} doesn't "
                        f'exist'
                    )
        return value
    return _callback


def safe_mkdir(dirpath: Path, descriptive_name: str = '') -> Path:
    """
    Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        sys.stderr.write(
            f'Path is empty: {descriptive_name if descriptive_name else ""}\n'
        )

    if dirpath.is_dir():
        return dirpath

    if dirpath.is_file():
        sys.stderr.write(f'{descriptive_name} {dirpath} is a file.\n')

    num_tries = 0
    max_tries = 10

    while not dirpath.exists():
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(str(dirpath))
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


def hash_sample_ids(sample_names: list[str]) -> str:
    """
    Return a unique hash string from a set of strings
    @param sample_names: set of strings
    @return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    h = hashlib.sha256(' '.join(sorted(sample_names)).encode()).hexdigest()[:38]
    return f'{h}_{len(sample_names)}'


def exists(path: Path | str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * cloud object
        * cloud URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    @param path: path to the file/directory/object/mt/ht
    @return: True if the object exists
    """
    path = cast(Path, to_path(path))

    # rstrip to ".mt/" -> ".mt"
    if any(str(path).rstrip('/').endswith(f'.{suf}') for suf in ['mt', 'ht']):
        path = path / '_SUCCESS'

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
