import logging
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd

import hail as hl
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, reference_path
from cpg_workflows.utils import can_reuse, exists
from gnomad.variant_qc.evaluation import compute_grouped_binned_ht
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg


def prepare_truth_sample_concordance(happy_vcf_path: str, high_confidence_only=True, reference_paths=None) -> hl.Table:
    """
    Prepare concordance data from a happy VCF file.

    :param happy_vcf_path: Path to the happy VCF file
    :param high_confidence_only: If True, filter to high confidence regions only
    :param reference_paths: Dict of dicts with structure:
                        {
                            'name': {
                                'path': 'gs://path/to/table.ht',
                                'keep': True  # True to keep regions, False to exclude
                            }
                        }
                        If None, uses default paths.

    :return: Hail Table with concordance annotations
    """
    # Load the concordance data
    happy_vcf = hl.import_vcf(happy_vcf_path, reference_genome='GRCh38', force_bgz=True)

    # Extract concordance metrics for the TRUTH and QUERY samples.
    fields = ['BD', 'BVT', 'BLT']
    samples = ['QUERY', 'TRUTH']

    annotations = {
        f'{field}_{sample}': hl.agg.filter(happy_vcf.s == sample, hl.agg.collect(happy_vcf[field]))[0]
        for field in fields
        for sample in samples
    }

    concordance = happy_vcf.annotate_rows(**annotations).rows()

    # Filter to high confidence regions if requested
    if high_confidence_only:
        # If no specific paths are provided, filter to high confidence NA12878 regions,
        # excluding segmental duplications and low complexity regions.
        if reference_paths is None:
            reference_paths = {
                'lcr': {
                    'path': reference_path('gnomad/lcr_intervals_ht'),
                    'keep': False,
                },
                'segdup': {
                    'path': reference_path('gnomad/seg_dup_intervals_ht'),
                    'keep': False,
                },
                'high_conf': {
                    'path': reference_path('na12878/regions_ht'),
                    'keep': True,
                },
            }

        # Filter to the specified regions.
        filter_expr = True
        for name, config in reference_paths.items():
            ref_table = hl.read_table(config['path'])
            if config['keep']:
                filter_expr = filter_expr & hl.is_defined(ref_table[concordance.locus])
            else:
                filter_expr = filter_expr & ~hl.is_defined(ref_table[concordance.locus])

        concordance = concordance.filter(filter_expr)

    return concordance


def create_binned_summary(
    vqsr_ht_path,
    happy_vcf_path: str,
    binned_summary_outpath: str,
    n_bins: int,
    fam_stats_ht_path: str | None,
    use_truth_sample_concordance: bool,
) -> hl.Table:
    """
    Create a binned summary of AS-VQSR variant quality scores with optional family and truth sample statistics.

    :param ht: Input Hail Table with AS-VQSR variant annotations
    :param happy_vcf_path: VCF path to pre-made table, required for use_truth_sample_concordance
    :param n_bins: Number of score bins to create
    :param fam_stats_ht: Optional family statistics table
    :param truth_sample_concordance: Optional truth sample concordance data
    :return: Binned summary table with aggregated statistics per score bin
    """
    # Ensure expected ht fields are present and appropriately named.

    # Load the VQSR ht
    vqsr_ht = hl.read_table(vqsr_ht_path)

    vqsr_ht = vqsr_ht.annotate(
        score=vqsr_ht.info.AS_VQSLOD,
        ac=vqsr_ht.info.AC,
        ac_raw=vqsr_ht.info.AC_raw,
        singleton=(vqsr_ht.info.AC_raw == 1),
        biallelic=(~vqsr_ht.was_split),
        positive_train_site=vqsr_ht.info.POSITIVE_TRAIN_SITE,
        negative_train_site=vqsr_ht.info.NEGATIVE_TRAIN_SITE,
        ac_qc_samples_unrelated_raw=vqsr_ht.info.AC_raw,
    )

    # Hacky addition to generate the truth_sample_concordance table within the function
    truth_sample_concordance = None

    if use_truth_sample_concordance:
        truth_sample_concordance = prepare_truth_sample_concordance(
            happy_vcf_path,
            high_confidence_only=True,
            reference_paths=None,
        )

    # If we have a truth sample, add those annotations too.
    if truth_sample_concordance is not None:
        vqsr_ht = vqsr_ht.annotate(
            TP=(truth_sample_concordance[vqsr_ht.key].BD_QUERY == "TP"),
            FP=(truth_sample_concordance[vqsr_ht.key].BD_QUERY == "FP"),
            FN=(truth_sample_concordance[vqsr_ht.key].BD_TRUTH == "FN"),
        )

    # Annotate with score bins, then group based on these.
    binned_ht = create_binned_ht(vqsr_ht, n_bins)
    grouped_binned_ht = compute_grouped_binned_ht(binned_ht)

    # Construct the aggregation functions; if no fam_stats_ht_path is provided
    # then exclude the de novo and transmitted singleton aggregations.

    # Hacky addition to read in fam_stats_ht if path provided
    if fam_stats_ht_path:
        fam_stats_ht = hl.read_table(fam_stats_ht_path)

    else:
        fam_stats_ht = None

    if fam_stats_ht is None:
        fam_stats_ht = hl.Table.parallelize(
            [],
            schema=hl.tstruct(
                locus=hl.tlocus(reference_genome='GRCh38'),
                alleles=hl.tarray(hl.tstr),
                n_de_novos_adj=hl.tint32,
                n_de_novos_raw=hl.tint32,
                n_transmitted_raw=hl.tint32,
                n_untransmitted_raw=hl.tint32,
                ac_parents_adj=hl.tint32,
                an_parents_adj=hl.tint32,
                ac_parents_raw=hl.tint32,
                an_parents_raw=hl.tint32,
            ),
        ).key_by('locus', 'alleles')

        all_aggregations = score_bin_agg(grouped_binned_ht, fam_stats_ht)
        aggregations = {k: v for k, v in all_aggregations.items() if 'de_novo' not in k and 'trans' not in k}
    else:
        aggregations = score_bin_agg(grouped_binned_ht, fam_stats_ht)

    # Add an aggregation for the number of biallelic variants.
    aggregations['n_biallelic'] = hl.agg.count_where(grouped_binned_ht._parent.biallelic)

    # If we have a truth sample, add concordance aggregations.
    if truth_sample_concordance is not None:
        aggregations['n_tp'] = hl.agg.count_where(grouped_binned_ht._parent.TP)
        aggregations['n_fp'] = hl.agg.count_where(grouped_binned_ht._parent.FP)
        aggregations['n_fn'] = hl.agg.count_where(grouped_binned_ht._parent.FN)

    # Calculate aggregated statistics for each score bin.
    agg_ht = grouped_binned_ht.aggregate(**aggregations)

    # For each unique bin_id (split into SNPs and indels), sum across the contigs.
    binned_summary_ht = agg_ht.group_by('bin_id', 'snv', 'bin').aggregate(
        min_score=hl.agg.min(agg_ht.min_score),
        max_score=hl.agg.max(agg_ht.max_score),
        **{
            x: hl.agg.sum(agg_ht[x])
            for x in agg_ht.row_value
            if x not in ['contig', 'bi_allelic', 'singleton', 'release_adj', 'min_score', 'max_score']
        },
    )
    binned_summary_ht.write(binned_summary_outpath, overwrite=True)

    return binned_summary_ht
