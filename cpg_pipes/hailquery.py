"""
Hail tables and matrix tables used as reference data.
"""
import asyncio
import functools
import logging
import operator
import os
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Union
import hail as hl

from .refdata import RefData
from .utils import safe_mkdir


logger = logging.getLogger(__file__)


def init_batch(billing_project: str, hail_bucket: Path | str):
    """
    Init Hail with Batch backend.
    """
    asyncio.get_event_loop().run_until_complete(
        hl.init_batch(
            default_reference=RefData.genome_build,
            billing_project=billing_project,
            remote_tmpdir=str(hail_bucket),
            token=os.environ['HAIL_TOKEN'],
        )
    )


def init_hail(name: str, local_tmp_dir: Path = None):
    """
    Initialise Hail, and set up a local directory for logs.
    @param name: name to prefix the log file
    @param local_tmp_dir: local directory to write Hail logs
    @return:
    """
    if not local_tmp_dir:
        local_tmp_dir = Path(tempfile.mkdtemp())

    timestamp = time.strftime('%Y%m%d-%H%M')
    hl_log = os.path.join(safe_mkdir(local_tmp_dir / 'log'), f'{name}-{timestamp}.log')
    hl.init(default_reference=RefData.genome_build, log=hl_log)
    return local_tmp_dir


def filter_low_conf_regions(
    mt: Union[hl.MatrixTable, hl.Table],
    refs: RefData,
    filter_lcr: bool = True,
    filter_segdup: bool = True,
    filter_telomeres_and_centromeres: bool = False,
    high_conf_regions: Optional[List[str]] = None,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter low-confidence regions.

    @param mt: MatrixTable or Table to filter
    @param filter_lcr: Whether to filter LCR regions
    @param filter_segdup: Whether to filter Segdup regions
    @param filter_telomeres_and_centromeres: Whether to filter telomeres and centromeres
    @param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    @return: MatrixTable or Table with low confidence regions removed
    """
    criteria = []
    if filter_lcr:
        lcr = hl.read_table(refs.lcr_intervals_ht)
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_segdup:
        segdup = hl.read_table(refs.seg_dup_intervals_ht)
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_telomeres_and_centromeres:
        telomeres_and_centromeres = hl.read_table(refs.tel_and_cent_ht)
        criteria.append(hl.is_missing(telomeres_and_centromeres[mt.locus]))

    if high_conf_regions is not None:
        for region in high_conf_regions:
            region = hl.import_locus_intervals(region)
            criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def get_truth_ht(refs: RefData) -> hl.Table:
    """
    Return a table with annotations from the latest version of the corresponding truth data.

    The following annotations are included:
        - hapmap
        - kgp_omni (1000 Genomes intersection Onni 2.5M array)
        - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
        - mills (Mills & Devine indels)

    @return: A table with the latest version of popular truth data annotations
    """
    return (
        hl.read_table(refs.hapmap_ht)
        .select(hapmap=True)
        .join(hl.read_table(refs.kgp_omni_ht).select(omni=True), how='outer')
        .join(hl.read_table(refs.kgp_hc_ht).select(kgp_phase1_hc=True), how='outer')
        .join(hl.read_table(refs.mills_ht).select(mills=True), how='outer')
        .repartition(200, shuffle=False)
        .persist()
    )
