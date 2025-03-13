import logging

import hail as hl

from cpg_utils import Path
from cpg_utils.config import config_retrieve, output_path
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import can_reuse
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX
from gnomad.resources.grch38.reference_data import (
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_adj,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    gen_anc_faf_max_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    region_flag_expr,
)
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict, make_freq_index_dict_from_meta


def run(
    vds_path: str,
    sample_qc_ht_path: str,
    relateds_to_drop_ht_path: str,
    infer_pop_ht_path: str,
    site_only_ht_path: str,
    out_ht_path: str,
):
    if can_reuse(out_ht_path):
        return

    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    vds = hl.vds.read_vds(str(vds_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))
    inferred_pop_ht = hl.read_table(str(infer_pop_ht_path))
    site_only_ht = hl.read_table(str(site_only_ht_path))

    logging.info('Generating frequency annotations...')
    freq_ht = frequency_annotations(
        vds,
        sample_qc_ht,
        relateds_to_drop_ht,
        inferred_pop_ht,
        site_only_ht,
        out_ht_path,
    )
    logging.info(f'Writing out frequency data to {out_ht_path}...')
    freq_ht.write(str(out_ht_path), overwrite=True)


def frequency_annotations(
    vds: hl.vds.VariantDataset,
    sample_qc_ht: hl.Table,
    relateds_to_drop_ht: hl.Table,
    inferred_pop_ht: hl.Table,
    site_only_ht: hl.Table,
    out_ht_path: str,
) -> hl.Table:
    """
    Generate frequency annotations (AF, AC, AN, InbreedingCoeff)
    """

    # Exome parsing
    if config_retrieve(['workflow'], ['sequencing_type']) == 'exome':
        # clear sample qc filters as this data has not been qc'd yet
        sample_qc_ht = sample_qc_ht.annotate(filters=hl.empty_set(hl.tstr))
        # sample_qc_ht has no females in it. This is causing errors in freq_meta_dict where we have no 'XY_adj'. Forcing
        # the one 'ambiguous' sex individual to be XY..
        sample_qc_ht = sample_qc_ht.annotate(
            X_karyotype=hl.if_else(
                (sample_qc_ht.X_karyotype == "ambiguous") & (sample_qc_ht.sex_karyotype == "ambiguous"),
                "X",
                sample_qc_ht.X_karyotype,
            ),
            Y_karyotype=hl.if_else(
                (sample_qc_ht.X_karyotype == "ambiguous") & (sample_qc_ht.sex_karyotype == "ambiguous"),
                "Y",
                sample_qc_ht.Y_karyotype,
            ),
            sex_karyotype=hl.if_else(
                (sample_qc_ht.X_karyotype == "ambiguous") & (sample_qc_ht.sex_karyotype == "ambiguous"),
                "XY",
                sample_qc_ht.sex_karyotype,
            ),
        )

    logging.info('Reading full sparse MT and metadata table...')

    logging.info('Splitting multiallelics')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logging.info('Densifying...')
    # The reason is that sparse matrix table has got NA records representing
    # the reference blocks, which affects the calculation of frequencies.
    # That's why we need to convert it to a "dense" representation, effectively
    # dropping reference blocks.
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    # Filter samples
    mt = mt.filter_cols(hl.len(sample_qc_ht[mt.col_key].filters) > 0, keep=False)
    mt = mt.filter_cols(hl.is_defined(relateds_to_drop_ht[mt.col_key]), keep=False)

    logging.info('Computing adj and sex adjusted genotypes...')
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, sample_qc_ht[mt.col_key].sex_karyotype),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    logging.info('Computing adj call rates...')
    mt_adj = mt.filter_entries(mt.adj)
    info_ht = mt_adj.annotate_rows(adj_gt_stats=hl.agg.call_stats(mt_adj.GT, mt_adj.alleles)).rows()

    logging.info('Generating frequency data...')
    mt = hl.variant_qc(mt)

    logging.info('Calculating InbreedingCoeff...')
    # NOTE: This is not the ideal location to calculate this, but added here
    # to avoid another densify.
    # The algorithm assumes all samples are unrelated:
    mt = mt.annotate_rows(InbreedingCoeff=hl.array([bi_allelic_site_inbreeding_expr(mt.GT)]))

    mt = annotate_labels(mt, inferred_pop_ht, sample_qc_ht)
    mt = _compute_filtering_af_and_popmax(mt)
    mt = mt.checkpoint(output_path('mt_faf_popmax.mt', category='tmp'), overwrite=True)
    # Currently have no Hail Tables with age data annotated on them, so unable to calculate age histograms
    # mt = _compute_age_hists(mt, sample_qc_ht)
    mt = mt.annotate_globals(freq_index_dict=make_freq_index_dict_from_meta(mt.freq_meta))

    # Annotate quality metrics histograms
    qual_hist_ht = _annotate_quality_metrics_hist(mt)
    mt = mt.annotate_rows(
        histograms=hl.struct(
            qual_hists=qual_hist_ht[mt.row_key].qual_hists,
            raw_qual_hists=qual_hist_ht[mt.row_key].raw_qual_hists,
        ),
    )

    freq_ht = mt.rows()
    freq_ht = freq_ht.annotate(info=site_only_ht[freq_ht.key].info)
    freq_ht = freq_ht.annotate(
        info=freq_ht.info.annotate(InbreedingCoeff=freq_ht.InbreedingCoeff),
    )
    freq_ht = freq_ht.annotate(
        region_flags=region_flag_expr(
            freq_ht,
            prob_regions={'lcr': lcr_intervals.ht(), 'segdup': seg_dup_intervals.ht()},
        ),
    )
    mono_allelic_flag_expr = (freq_ht.freq[0].AC > 0) & (freq_ht.freq[1].AF == 1)
    freq_ht = freq_ht.annotate(info=freq_ht.info.annotate(monoallelic=mono_allelic_flag_expr))
    freq_ht = freq_ht.annotate(**freq_ht.variant_qc)
    # freq_ht = freq_ht.drop('variant_qc')
    freq_ht = freq_ht.annotate(**info_ht[freq_ht.locus, freq_ht.alleles].select('adj_gt_stats'))

    return freq_ht


def annotate_labels(mt: hl.MatrixTable, inferred_pop_ht: hl.Table, sample_qc_ht: hl.Table) -> hl.MatrixTable:
    # prepare_gnomad_v4_variants_helper requires ancestry to be annotated
    print(f'{mt.s.collect()}')
    print(f'{inferred_pop_ht.s.collect()}')
    print(f'{sample_qc_ht.s.collect()}')
    print(mt.show())
    print(inferred_pop_ht.show())
    print(sample_qc_ht.show())
    mt = mt.annotate_cols(gen_anc=inferred_pop_ht[mt.s].pop)
    # prepare_gnomad_v4_variants_helper requires sex to be annotated
    mt = mt.annotate_cols(sex=sample_qc_ht[mt.s].sex_karyotype)
    # mt = mt.annotate_cols(subset='tenk10k')
    mt = annotate_freq(
        mt,
        sex_expr=mt.sex,
        additional_strata_expr=[{'gen_anc': mt.gen_anc}],
        pop_expr=mt.gen_anc,
    )
    return mt


def _compute_age_hists(mt: hl.MatrixTable, sample_qc_ht: hl.Table) -> hl.MatrixTable:
    logging.info('Computing age histograms for each variant...')
    try:
        mt = mt.annotate_cols(age=hl.float64(sample_qc_ht[mt.col_key].age))
    except AttributeError:
        pass
    else:
        mt = mt.annotate_rows(
            **age_hists_expr(
                mt.adj,
                mt.GT,
                mt.age,
            ),
        )
        # Compute callset-wide age histogram global
        mt = mt.annotate_globals(
            age_distribution=mt.aggregate_cols(
                hl.agg.hist(
                    mt.age,
                    30,
                    80,
                    10,
                ),
            ),
        )
    return mt


def _compute_filtering_af_and_popmax(mt: hl.MatrixTable) -> hl.MatrixTable:
    logging.info('Computing filtering allele frequencies and popmax...')
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX, pop_label='gen_anc')
    mt = mt.select_rows(
        'InbreedingCoeff',
        'freq',
        'rsid',
        'variant_qc',
        faf=faf,
        popmax=pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX),
    )
    mt = mt.annotate_globals(
        faf_meta=faf_meta,
        # exome data has 'Other' while genome data has 'oth'
        faf_index_dict=make_faf_index_dict(faf_meta, pops=['Europe', 'oth', 'Other']),
    )
    mt = mt.annotate_rows(
        popmax=mt.popmax.annotate(
            faf95=mt.faf[mt.faf_meta.index(lambda x: x.values() == ['adj', mt.popmax.pop])].faf95,
        ),
    )
    mt = mt.annotate_rows(fafmax=gen_anc_faf_max_expr(faf=mt.faf, faf_meta=mt.faf_meta, pop_label='gen_anc'))

    # Populating 'filters' field with empty set for now
    mt = mt.annotate_rows(filters=hl.empty_set(hl.tstr))
    return mt


def _annotate_quality_metrics_hist(mt: hl.MatrixTable) -> hl.Table:
    logging.info('Annotating quality metrics histograms...')
    # NOTE: these are performed here as the quality metrics histograms
    # also require densifying
    mt = mt.annotate_rows(qual_hists=qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj))
    ht = mt.rows()
    ht = ht.annotate(
        qual_hists=hl.Struct(**{i.replace('_adj', ''): ht.qual_hists[i] for i in ht.qual_hists if '_adj' in i}),
        raw_qual_hists=hl.Struct(**{i: ht.qual_hists[i] for i in ht.qual_hists if '_adj' not in i}),
    )
    return ht
