import logging

import hail as hl

from cpg_utils import Path
from cpg_utils.config import config_retrieve, output_path
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import can_reuse
from gnomad.resources.grch38.reference_data import (
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    gen_anc_faf_max_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    region_flag_expr,
)
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict_from_meta

POPS_TO_REMOVE_FOR_POPMAX = {'Unclassified'}


def run(
    vds_path: str,
    sample_qc_ht_path: str,
    relateds_to_drop_ht_path: str,
    vqsr_ht_path: str,
    out_ht_path: str,
):
    if can_reuse(out_ht_path):
        return

    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    vds = hl.vds.read_vds(str(vds_path))
    sample_qc_ht = hl.read_table(str(sample_qc_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))
    vqsr_ht = hl.read_table(str(vqsr_ht_path))

    logging.info('Generating frequency annotations...')
    freq_ht = frequency_annotations(vds, sample_qc_ht, relateds_to_drop_ht, vqsr_ht)
    logging.info(f'Writing out frequency data to {out_ht_path}...')
    freq_ht.write(str(out_ht_path), overwrite=True)


def frequency_annotations(
    vds: hl.vds.VariantDataset,
    sample_qc_ht: hl.Table,
    relateds_to_drop_ht: hl.Table,
    vqsr_ht: hl.Table,
) -> hl.Table:
    """
    Generate frequency annotations (AF, AC, AN, InbreedingCoeff)
    """

    # Prepare the dense matrix table.
    logging.info('Splitting multiallelics...')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logging.info('Densifying...')
    mt = hl.vds.to_dense_mt(vds)

    # Filter and annotate samples.
    # Note: removing samples here may leave an otherwise dense MT with hom-ref genotypes
    logging.info('Removing filtered samples...')
    mt = mt.filter_cols(hl.len(sample_qc_ht[mt.col_key].filters) > 0, keep=False)
    mt = mt.filter_cols(hl.is_defined(relateds_to_drop_ht[mt.col_key]), keep=False)

    logging.info('Annotating sample metadata...')
    mt = mt.annotate_cols(
        gen_anc=sample_qc_ht[mt.s].population,
        sex=sample_qc_ht[mt.s].sex_karyotype,
        subset=sample_qc_ht[mt.s].subset if "subset" in list(sample_qc_ht.row) else hl.missing('str'),
    )

    logging.info('Computing quality and sex adjusted genotypes...')
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, sample_qc_ht[mt.col_key].sex_karyotype),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    mt = mt.checkpoint(output_path('filtered_checkpoint.mt', category='tmp'), overwrite=True)

    # Annotate variant QC metrics.
    logging.info('Annotating variant QC metrics...')
    logging.info('Base variant QC...')
    mt = hl.variant_qc(mt)

    logging.info('Allele-specific statistics...')
    mt = mt.annotate_rows(info=vqsr_ht[mt.row_key].info)
    mt = mt.annotate_rows(AS_lowqual=vqsr_ht[mt.row_key].AS_lowqual)

    logging.info('Inbreeding coefficient...')
    mt = mt.annotate_rows(inbreeding_coeff=hl.array([bi_allelic_site_inbreeding_expr(mt.GT)]))

    logging.info('VQSR filters...')
    mt = mt.annotate_rows(site_vqsr_filters=vqsr_ht[mt.row_key].filters)
    mt = mt.annotate_rows(as_vqsr_filters=vqsr_ht[mt.row_key].info.AS_FilterStatus)

    logging.info('Region flags...')
    mt = mt.annotate_rows(
        region_flags=region_flag_expr(
            mt,
            prob_regions={'lcr': lcr_intervals.ht(), 'segdup': seg_dup_intervals.ht()},
        ),
    )

    mt = mt.checkpoint(output_path('variant_qc_checkpoint.mt', category='tmp'), overwrite=True)

    # Compute allele frequencies, stratified by genetic ancestry, sex, and data subset (if applicable).
    logging.info('Computing stratified allele frequencies...')
    mt = annotate_freq(
        mt,
        sex_expr=mt.sex,
        pop_expr=mt.gen_anc,
        additional_strata_expr=[
            {'subset': mt.subset},
            {'subset': mt.subset, 'pop': mt.gen_anc},
            {'subset': mt.subset, 'pop': mt.gen_anc, 'sex': mt.sex},
        ],
    )
    mt = mt.checkpoint(output_path('stratified_frequencies.mt', category='tmp'), overwrite=True)
    mt = mt.annotate_globals(freq_index_dict=make_freq_index_dict_from_meta(mt.freq_meta))
    mt = _compute_filtering_af_and_popmax(mt)
    mt = mt.annotate_rows(monoallelic=(mt.freq[0].AC > 0) & (mt.freq[1].AF == 1))

    # Compute variant histograms.
    logging.info('Computing variant histograms...')
    mt = _compute_age_hists(mt, sample_qc_ht)
    qual_hist_ht = _annotate_quality_metrics_hist(mt)
    mt = mt.annotate_rows(
        histograms=hl.struct(
            qual_hists=qual_hist_ht[mt.row_key].qual_hists,
            raw_qual_hists=qual_hist_ht[mt.row_key].raw_qual_hists,
        ),
    )

    return mt.rows()


def _compute_age_hists(mt: hl.MatrixTable, sample_qc_ht: hl.Table) -> hl.MatrixTable:
    logging.info('Computing carrier age histograms for each variant...')
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
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX, pop_label='pop')
    popmax = pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX, pop_label='pop')
    mt = mt.annotate_rows(
        faf=faf,
        popmax=popmax,
    )

    mt = mt.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(
            faf_meta,
            pops=list(mt.aggregate_cols(hl.agg.collect_as_set(mt.gen_anc))),  # unique pop labels,
        ),
    )

    # TO DO: adjust fafmax expression to correctly handle the case when AF != MAF.
    mt = mt.annotate_rows(
        fafmax=gen_anc_faf_max_expr(faf=mt.faf, faf_meta=mt.faf_meta, pop_label='pop'),
        popmax=mt.popmax.annotate(
            faf95=mt.faf[
                mt.faf_meta.index(
                    lambda x: (
                        (x.key_set() == hl.set(['group', 'pop'])) & (x['group'] == 'adj') & (x['pop'] == mt.popmax.pop)
                    ),
                )
            ].faf95,
        ),
    )

    return mt


def _annotate_quality_metrics_hist(mt: hl.MatrixTable) -> hl.Table:
    logging.info('Computing variant quality metrics histograms...')
    # NOTE: these are performed here as the quality metrics histograms
    # also require densifying
    mt = mt.annotate_rows(qual_hists=qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj))
    ht = mt.rows()
    ht = ht.annotate(
        qual_hists=hl.Struct(**{i.replace('_adj', ''): ht.qual_hists[i] for i in ht.qual_hists if '_adj' in i}),
        raw_qual_hists=hl.Struct(**{i: ht.qual_hists[i] for i in ht.qual_hists if '_adj' not in i}),
    )
    return ht
