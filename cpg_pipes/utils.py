"""
Various utility functions and constants.
"""

import hashlib
import os
import sys
import time
from dataclasses import dataclass
from enum import Enum
from os.path import isdir, isfile, exists
from typing import Any, Callable, Dict, Optional, Iterable, Set

import click

from cpg_pipes import __name__ as package_name

# Default reference genome build.
from cpg_pipes.buckets import file_exists
DEFAULT_REF = 'GRCh38'

# Packages to install on a dataproc cluster, to use with the dataproc wrapper.
DATAPROC_PACKAGES = [
    'cpg-pipes',
    'cpg-gnomad',   # github.com/populationgenomics/gnomad_methods
    'seqr-loader',  # github.com/populationgenomics/hail-elasticsearch-pipelines
    'elasticsearch',
    'click',
    'google',
    'fsspec',
    'sklearn',
    'gcloud',
    'selenium',
]

# Location of python scripts to be called directly from command line.
SCRIPTS_DIR = 'scripts'

# Location of Hail Query scripts, to use with the dataproc wrapper.
QUERY_SCRIPTS_DIR = 'query_scripts'

# This python package name.
PACKAGE_DIR = package_name


class Namespace(Enum):
    """
    CPG storage namespace. See for more details on storage policies:
    https://github.com/populationgenomics/team-docs/tree/main/storage_policies#main-vs-test
    """
    MAIN = 'main'  # production runs: read from main, write to main
    TEST = 'test'  # runs that make test data: read from test, write to test
    TMP = 'tmp'    # read from test, write to tmp


class AnalysisType(Enum):
    """
    Types of sample-metadata analysis. Corresponds to Analysis types:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L4-L11

    Re-defined in utils to allow using the class in type hints without 
    importing sample-metadata.
    """
    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'


@dataclass
class Analysis:
    """
    Analysis SampleMetadata DB entry. 

    See sample-metadata for more details: https://github.com/populationgenomics/sample-metadata

    Defined in utils to allow using the class in type hints without importing 
    sample-metadata.
    """
    id: int
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]


class Sequence:
    """
    Sequence SampleMetadata DB entry.

    See sample-metadata for more details: https://github.com/populationgenomics/sample-metadata

    Defined in utils to allow using the class in type hints without importing 
    sample-metadata.
    """
    def __init__(self, id, sample_id, meta, smdb):
        self.id = id
        self.sample_id = sample_id
        self.meta = meta
        self.smdb = smdb
    
    @staticmethod
    def parse(data: Dict, smdb):
        return Sequence(data['id'], data['sample_id'], data['meta'], smdb)


def get_validation_callback(
    ext: str = None,
    must_exist: bool = False,
    accompanying_metadata_suffix: str = None,
) -> Callable:
    """
    Get callback for Click parameters validation
    :param ext: check that the path has the expected extension
    :param must_exist: check that the input file/object/directory exists
    :param accompanying_metadata_suffix: checks that a file at the same location but
    with a different suffix also exists (e.g. genomes.mt and genomes.metadata.ht)
    :return: a callback suitable for Click parameter initialization
    """
    def callback(_: click.Context, param: click.Option, value: Any):
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
            if not file_exists(value):
                raise click.BadParameter(f"{value} doesn't exist or incomplete")
            if accompanying_metadata_suffix:
                accompanying_metadata_fpath = (
                    os.path.splitext(value)[0] + accompanying_metadata_suffix
                )
                if not file_exists(accompanying_metadata_fpath):
                    raise click.BadParameter(
                        f"An accompanying file {accompanying_metadata_fpath} doesn't "
                        f'exist'
                    )
        return value

    return callback


def safe_mkdir(dirpath: str, descriptive_name: str = '') -> str:
    """
    Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        sys.stderr.write(
            f'Path is empty: {descriptive_name if descriptive_name else ""}\n'
        )

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        sys.stderr.write(descriptive_name + ' ' + dirpath + ' is a file.\n')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


def hash_sample_ids(sample_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a set of strings
    :param sample_names: set of strings
    :return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    return hashlib.sha256(' '.join(sorted(sample_names)).encode()).hexdigest()[:32]
