"""
Utility functions for the joint_calling module.
"""

import os
import subprocess
import tempfile
import logging
import sys
import time
import hashlib
from dataclasses import dataclass
from enum import Enum
from os.path import isdir, isfile, exists, join, basename
from typing import Any, Callable, Dict, Optional, Union, Iterable, Set
import yaml
import pandas as pd
import hail as hl
import click
from google.cloud import storage
from cpg_pipes import _version, get_package_path
from cpg_pipes import __name__ as package_name


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


# Default reference build.
DEFAULT_REF = 'GRCh38'

# Packages to install on a dataproc cluster, to use with the dataproc wrapper.
DATAPROC_PACKAGES = [
    'cpg-pipes',
    'cpg-gnomad',
    'seqr-loader',
    'elasticsearch',
    'click',
    'google',
    'slackclient',
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


def init_hail(name: str, local_tmp_dir: str = None):
    """
    Initialise Hail, and set up a local directory for logs.
    :param name: name to prefix the log file
    :param local_tmp_dir: local directory to write Hail logs
    :return:
    """
    if not local_tmp_dir:
        local_tmp_dir = tempfile.mkdtemp()

    timestamp = time.strftime('%Y%m%d-%H%M')
    hl_log = os.path.join(
        safe_mkdir(os.path.join(local_tmp_dir, 'log')), f'{name}-{timestamp}.log'
    )
    hl.init(default_reference=DEFAULT_REF, log=hl_log)
    logger.info(f'Running joint-calling version {_version.__version__}')
    return local_tmp_dir


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


def gs_cache_file(fpath: str, local_tmp_dir: str) -> str:
    """
    :param fpath: local or a `gs://` path. If the latter, the file
        will be downloaded and cached if local_tmp_dir is provided,
        the local path will be returned
    :param local_tmp_dir: a local directory to cache files downloaded
        from Google Storage
    :return: file path
    """
    if fpath.startswith('gs://'):
        fname = (
            os.path.basename(fpath) + '_' + hashlib.md5(fpath.encode()).hexdigest()[:6]
        )
        local_fpath = os.path.join(local_tmp_dir, fname)
        if not exists(local_fpath):
            bucket = fpath.replace('gs://', '').split('/')[0]
            path = fpath.replace('gs://', '').split('/', maxsplit=1)[1]
            gs = storage.Client()
            blob = gs.get_bucket(bucket).get_blob(path)
            if blob:
                blob.download_to_filename(local_fpath)
    else:
        local_fpath = fpath
    return local_fpath


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


def get_mt(
    mt_path: str,
    split: bool = False,
    hard_filtered_samples_to_remove_ht: hl.Table = None,
    meta_ht: hl.Table = None,
    add_meta: bool = False,
    release_only: bool = False,
    passing_sites_only: bool = False,
) -> hl.MatrixTable:
    """
    Wrapper function to get data with desired filtering and metadata annotations
    :param mt_path: path to the MatrixTable
    :param split:
        Split multiallelics and convert local-allele LGT/LA fields to GT.
        Note: this will perform a split on the MT rather than grab an already split MT
    :param hard_filtered_samples_to_remove_ht:
        table with samples to remove
        (only relevant after sample QC that produces a table with samples failed
        filtering)
    :param meta_ht: table with meta-information generated by sample QC
    :param add_meta: whether to add metadata to MT in 'meta' column
    :param release_only: whether to filter the MT to only samples available for
        release (can only be used if metadata is present)
    :param passing_sites_only: whether to filter the MT to only variants with
        nothing in the filter field (e.g. passing soft filters)
    :return: MatrixTable with chosen annotations and filters
    """
    mt = hl.read_matrix_table(mt_path)

    if passing_sites_only:
        try:
            mt = mt.filter_rows(mt.filters.length() == 0)
        except AttributeError:
            pass

    if hard_filtered_samples_to_remove_ht is not None:
        mt = mt.filter_cols(
            hl.is_missing(hard_filtered_samples_to_remove_ht[mt.col_key])
        )

    if add_meta:
        assert meta_ht is not None
        mt = mt.annotate_cols(meta=meta_ht[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

    elif release_only:
        assert meta_ht is not None
        mt = mt.filter_cols(meta_ht[mt.col_key].release)

    if split:
        mt = mt.annotate_rows(
            n_unsplit_alleles=hl.len(mt.alleles),
            mixed_site=(hl.len(mt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]),
        )
        # Will use GT instead of LGT
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    return mt


def get_vqsr_filters_path(
    work_bucket: str,
    model_id: str,
    split: bool = True,
    finalized: bool = False,
) -> str:
    """
    Gets the specified VQSR filtering annotation resource.
    :param work_bucket: bucket
    :param model_id: VQSR filtering model id
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file path
    """
    return join(
        work_bucket,
        f'filtering/{model_id}'
        f'{".finalized" if finalized else ""}'
        f'{".split" if split else ""}'
        f'.ht',
    )


def get_filter_cutoffs(
    provided_filter_cutoffs_path: Optional[str] = None,
) -> Dict:
    """
    :provided_filter_cutoffs_path: optional, a path to a YAML file with cutoffs.
    Can sit on a bucket. If not provided, a default one from the package will be used.
    gets the a default one within the package
    :return: a Dict with cutoffs
    """
    if provided_filter_cutoffs_path:
        assert file_exists(provided_filter_cutoffs_path), provided_filter_cutoffs_path
        path = provided_filter_cutoffs_path
    else:
        path = join(get_package_path(), 'filter_cutoffs.yaml')

    if path.startswith('gs://'):
        contents = subprocess.check_output(['gsutil', 'cat', path])
        filter_cutoffs_d = yaml.safe_load(contents)
    else:
        with open(path) as f:
            filter_cutoffs_d = yaml.safe_load(f)

    return filter_cutoffs_d


def parse_input_metadata(
    meta_csv_path: str,
    local_tmp_dir: str,
    out_ht_path: Optional[str] = None,
) -> hl.Table:
    """
    Parse KCCG metadata (continental_pop and picard metrics)
    """
    local_csv_path = join(local_tmp_dir, basename(meta_csv_path))
    gsutil_cp(meta_csv_path, local_csv_path)
    df = pd.read_table(local_csv_path)
    ht = hl.Table.from_pandas(df).key_by('s')
    if out_ht_path:
        ht = ht.checkpoint(out_ht_path, overwrite=True)
    return ht


def hash_sample_ids(sample_names: Iterable[str]) -> str:
    """
    Return a unique hash string from a set of strings
    :param sample_names: set of strings
    :return: a string hash
    """
    for sn in sample_names:
        assert ' ' not in sn, sn
    return hashlib.sha256(' '.join(sorted(sample_names)).encode()).hexdigest()[:32]


def gsutil_cp(
    src_path: str,
    dst_path: str,
    disable_check_hashes: bool = False,
    recursive: bool = False,
    quiet: bool = False,
):
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
