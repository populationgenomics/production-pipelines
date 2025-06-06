"""
Utility functions and constants.
"""

import logging
import re
import string
import sys
import time
import traceback
import unicodedata
from functools import lru_cache
from itertools import chain, islice
from os.path import basename, dirname, join
from random import choices
from typing import Union, cast

import coloredlogs

import hail as hl
from hailtop.batch import ResourceFile, ResourceGroup

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import get_batch

DEFAULT_LOG_FORMAT = config_retrieve(
    ['workflow', 'logger', 'default_format'],
    '%(asctime)s - %(name)s - %(pathname)s: %(lineno)d - %(levelname)s - %(message)s',
)
LOGGERS: dict[str, logging.Logger] = {}

VCF_BGZ = 'vcf.bgz'
VCF_BGZ_TBI = 'vcf.bgz.tbi'
VCF_GZ = 'vcf.gz'
VCF_GZ_TBI = 'vcf.gz.tbi'


def get_logger(
    logger_name: str = 'cpg_workflows',
    log_level: int = logging.INFO,
    fmt_string: str = DEFAULT_LOG_FORMAT,
) -> logging.Logger:
    """
    creates a logger instance (so as not to use the root logger)
    Args:
        logger_name (str):
        log_level (int): logging level, defaults to INFO. Can be overridden by config
        fmt_string (str): format string for this logger, defaults to DEFAULT_LOG_FORMAT
    Returns:
        a logger instance, if required create it first
    """

    if logger_name not in LOGGERS:

        # allow a log-level & format override on a name basis
        log_level = config_retrieve(['workflow', 'logger', logger_name, 'level'], log_level)
        fmt_string = config_retrieve(['workflow', 'logger', logger_name, 'format'], fmt_string)

        # create a named logger
        new_logger = logging.getLogger(logger_name)
        new_logger.setLevel(log_level)

        # unless otherwise specified, use coloredlogs
        if config_retrieve(['workflow', 'logger', logger_name, 'use_colored_logs'], False):
            coloredlogs.install(level=log_level, fmt=fmt_string, logger=new_logger)

        # create a stream handler to write output
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(log_level)

        # create format string for messages
        formatter = logging.Formatter(fmt_string)
        stream_handler.setFormatter(formatter)

        # set the logger to use this handler
        new_logger.addHandler(stream_handler)

        LOGGERS[logger_name] = new_logger

    return LOGGERS[logger_name]


def chunks(iterable, chunk_size):
    """
    Yield successive n-sized chunks from an iterable

    Args:
        iterable (): any iterable - tuple, str, list, set
        chunk_size (): size of intervals to return

    Returns:
        intervals of requested size across the collection
    """

    if isinstance(iterable, set):
        iterable = list(iterable)

    for i in range(0, len(iterable), chunk_size):
        yield iterable[i : (i + chunk_size)]


def generator_chunks(generator, size):
    """
    Iterates across a generator, returning specifically sized chunks

    Args:
        generator (): any generator or method implementing yield
        size (): size of iterator to return

    Returns:
        a subset of the generator results
    """
    iterator = iter(generator)
    for first in iterator:
        yield list(chain([first], islice(iterator, size - 1)))


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


def checkpoint_hail(
    t: hl.Table | hl.MatrixTable,
    file_name: str,
    checkpoint_prefix: str | None = None,
    allow_reuse=False,
):
    """
    checkpoint method
    provide with a path and a prefix (GCP directory, can be None)
    allow_reuse sets whether the checkpoint can be reused - we
    typically want to avoid reuse, as it means we're continuing a previous
    failure from an unknown state

    Args:
        t (hl.Table | hl.MatrixTable):
        file_name (str): name for this checkpoint
        checkpoint_prefix (str): path to the checkpoint directory
        allow_reuse (bool): whether to permit reuse of an existing checkpoint
    """

    # drop the schema here
    t.describe()

    # log the current number of partitions
    logging.info(f'Checkpointing object as {t.n_partitions()} partitions')

    if checkpoint_prefix is None:
        return t

    path = join(checkpoint_prefix, file_name)
    if can_reuse(path) and allow_reuse:
        logging.info(f'Re-using {path}')
        return read_hail(path)

    logging.info(f'Checkpointing {path}')
    return t.checkpoint(path, overwrite=True)


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
            res = check_exists_path(path)

        # a failure to detect the parent folder causes a crash
        # instead stick to a core responsibility -
        # existence = False
        except FileNotFoundError as fnfe:
            logging.error(f'Failed checking {path}')
            logging.error(f'{fnfe}')
            return False
        except BaseException:
            traceback.print_exc()
            logging.error(f'Failed checking {path}')
            sys.exit(1)
        logging.debug(f'Checked {path} [' + ('exists' if res else 'missing') + ']')
        return res

    return check_exists_path(path)


def check_exists_path(test_path: Path) -> bool:
    """
    Check whether a path exists using a cached per-directory listing.
    NB. reversion to Strings prevents a get call, which is typically
    forbidden to local users - this prevents this method being used in the
    metamist audit processes
    """
    return basename(str(test_path)) in get_contents_of_path(dirname(str(test_path)))


@lru_cache
def get_contents_of_path(test_path: str) -> set[str]:
    """
    Get the contents of a GCS path, returning non-complete paths, eg:

        get_contents_of_path('gs://my-bucket/my-dir/')
        'my-file.txt'

    """
    return {f.name for f in to_path(test_path.rstrip('/')).iterdir()}


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
        rand_bit = ''.join(choices(string.ascii_uppercase + string.digits, k=rand_suffix_len))
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

    @param rich_id_map: map used to replace sequencing groups, e.g. {'CPGAA': 'CPGAA|EXTID'}
    @param file_names: file names and Hail Batch Resource files where to replace IDs
    @return: bash command that does replacement
    """
    cmd = ''
    for sgid, rich_sgid in rich_id_map.items():
        for fname in file_names:
            cmd += f'sed -iBAK \'s/{sgid}/{rich_sgid}/g\' {fname}'
            cmd += '\n'
    return cmd


def tshirt_mt_sizing(sequencing_type: str, cohort_size: int) -> int:
    """
    Some way of taking the details we have (#SGs, sequencing type)
    and producing an estimate (with padding) of the MT size on disc
    used to determine VM provision during ES export and Talos

    Args:
        sequencing_type ():
        cohort_size ():

    Returns:
        str, the value for job.storage(X)
    """

    # allow for an override from config
    if preset := config_retrieve(['workflow', 'es_storage'], False):
        return preset

    if (sequencing_type == 'genome' and cohort_size < 100) or (sequencing_type == 'exome' and cohort_size < 1000):
        return 50
    return 500


def get_intervals_from_bed(intervals_path: Path) -> list[str]:
    """
    Read genomic intervals from a bed file.
    Increment the start position of each interval by 1 to match the 1-based
    coordinate system used by GATK.

    Returns a list of interval strings in the format 'chrN:start-end'.
    """
    with intervals_path.open('r') as f:
        intervals = []
        for line in f:
            chrom, start, end = line.strip().split('\t')
            intervals.append(f'{chrom}:{int(start)+1}-{end}')
    return intervals


@lru_cache(2)
def get_all_fragments_from_manifest(manifest_file: Path) -> list[ResourceGroup]:
    """
    read the manifest file, and return all the fragment resources as an ordered list
    this is a cached method as we don't want to localise every fragment once per task

    Args:
        manifest_file ():

    Returns:
        an ordered list of all the fragment VCFs and corresponding indices
    """

    resource_objects: list[ResourceGroup] = []
    manifest_folder: Path = manifest_file.parent
    with manifest_file.open() as f:
        for line in f:
            vcf_path = manifest_folder / line.strip()
            resource_objects.append(
                get_batch().read_input_group(**{VCF_GZ: vcf_path, VCF_GZ_TBI: f'{vcf_path}.tbi'}),
            )
    return resource_objects


ExpectedResultT = Union[Path, dict[str, Path], dict[str, str], dict[str, Path | str], str, None]
