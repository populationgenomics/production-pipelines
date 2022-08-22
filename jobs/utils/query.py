"""
Hail tables and matrix tables used as reference data.
"""
import functools
import logging
import operator
from typing import List, Optional, Union

import hail as hl

from cpg_utils.hail_batch import reference_path


def filter_low_conf_regions(
    mt: Union[hl.MatrixTable, hl.Table],
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
    @param filter_telomeres_and_centromeres: Whether to filter telomeres and
        centromeres
    @param high_conf_regions: Paths to set of high confidence regions to restrict to
        (union of regions)
    @return: MatrixTable or Table with low confidence regions removed
    """
    criteria = []
    if filter_lcr:
        lcr = hl.read_table(str(reference_path('gnomad/lcr_intervals_ht')))
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_segdup:
        segdup = hl.read_table(str(reference_path('gnomad/seg_dup_intervals_ht')))
        criteria.append(hl.is_missing(segdup[mt.locus]))

    if filter_telomeres_and_centromeres:
        telomeres_and_centromeres = hl.read_table(
            str(reference_path('gnomad/tel_and_cent_ht'))
        )
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


def get_truth_ht() -> hl.Table:
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
        hl.read_table(str(reference_path('gnomad/hapmap_ht')))
        .select(hapmap=True)
        .join(
            hl.read_table(str(reference_path('gnomad/kgp_omni_ht'))).select(omni=True),
            how='outer',
        )
        .join(
            hl.read_table(str(reference_path('gnomad/kgp_hc_ht'))).select(
                kgp_phase1_hc=True
            ),
            how='outer',
        )
        .join(
            hl.read_table(str(reference_path('gnomad/mills_ht'))).select(mills=True),
            how='outer',
        )
        .repartition(200, shuffle=False)
        .persist()
    )
