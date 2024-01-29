"""
Impute sex.
Add soft filters for samples.
"""

import logging

import hail as hl
from hail.vds.variant_dataset import VariantDataset
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build, reference_path
from gnomad.sample_qc.pipeline import annotate_sex

from cpg_workflows.inputs import get_cohort
from cpg_workflows.utils import can_reuse


def run(
    vds_path: Path,
    out_sample_qc_ht_path: Path,
    tmp_prefix: Path,
):
    if can_reuse(out_sample_qc_ht_path, overwrite=True):
        return hl.read_table(str(out_sample_qc_ht_path))

    ht = initialise_sample_table()

    vds = hl.vds.read_vds(str(vds_path))

    # Remove centromeres and telomeres:
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    if tel_cent_ht.count() > 0:
        vds = hl.vds.filter_intervals(vds, tel_cent_ht, keep=False)

    # Run Hail sample-QC stats:
    sqc_ht_path = tmp_prefix / 'sample_qc.ht'
    if can_reuse(sqc_ht_path, overwrite=True):
        sqc_ht = hl.read_table(str(sqc_ht_path))
    else:
        # Filter to autosomes:
        autosome_vds = hl.vds.filter_chromosomes(
            vds, keep=[f'chr{chrom}' for chrom in range(1, 23)]
        )
        sqc_ht = hl.vds.sample_qc(autosome_vds)
        sqc_ht = sqc_ht.checkpoint(str(sqc_ht_path), overwrite=True)
    ht = ht.annotate(sample_qc=sqc_ht[ht.s])
    logging.info('Sample QC table:')
    ht.describe()

    logging.info('Run sex imputation')
    sex_ht = impute_sex(vds, ht, tmp_prefix)
    ht = ht.annotate(**sex_ht[ht.s])

    logging.info('Adding soft filters')
    ht = add_soft_filters(ht)
    ht.checkpoint(str(out_sample_qc_ht_path), overwrite=True)


def initialise_sample_table() -> hl.Table:
    """
    Export the cohort into a sample-level Hail Table.
    """
    training_pop = get_config()['large_cohort'].get('training_pop')
    entries = [
        {
            's': s.id,
            'external_id': s.external_id,
            'dataset': s.dataset.name,
            'gvcf': str(s.gvcf.path),
            'sex': s.pedigree.sex.value,
            'training_pop': s.meta.get(training_pop) if training_pop else None,
            'superpopulation': s.meta.get('Superpopulation name'),
            'population': s.meta.get('Population name'),
        }
        for s in get_cohort().get_sequencing_groups()
        if s.gvcf
    ]
    if not entries:
        raise ValueError('No samples with GVCFs found')
    t = 'array<struct{s: str, external_id: str, dataset: str, gvcf: str, sex: int, training_pop: str, superpopulation: str, population: str}>'
    ht = hl.Table.parallelize(hl.literal(entries, t), key='s')
    return ht


def generate_sex_coverage_mt(
    vds: hl.vds.VariantDataset,
    calling_intervals: hl.Table,
    tmp_prefix: Path,
):
    """
    generate the pre-computed MT of coverage intervals, to supply to annotate_sex
    Args:
        vds ():
        calling_intervals ():
        tmp_prefix ():

    Returns:

    """
    key_dtype = calling_intervals.key.dtype
    if (
            len(key_dtype) != 1 or
            not isinstance(calling_intervals.key[0].dtype, hl.tinterval) or
            calling_intervals.key[0].dtype.point_type != vds.reference_data.locus.dtype
    ):
        raise ValueError(f'"impute_sex_chromosome_ploidy": expect calling_intervals to be list of intervals or '
                         f'table with single key of type interval<locus>, found table with key: {key_dtype}')

    rg = vds.reference_data.locus.dtype.reference_genome
    par_boundaries = []
    for par_interval in rg.par:
        par_boundaries.append(par_interval.start)
        par_boundaries.append(par_interval.end)

    # segment on PAR interval boundaries
    calling_intervals = hl.segment_intervals(calling_intervals, par_boundaries)

    # remove intervals overlapping PAR
    calling_intervals = calling_intervals.filter(
        hl.all(lambda x: ~x.overlaps(calling_intervals.interval), hl.literal(rg.par)))

    # checkpoint for efficient multiple downstream usages
    checkpoint_path = tmp_prefix / 'sample_qc' / 'calling_intervals_partial.ht'
    if can_reuse(checkpoint_path):
        calling_intervals = hl.read_matrix_table(str(checkpoint_path))
    else:
        calling_intervals = calling_intervals.checkpoint(str(checkpoint_path), overwrite=True)

    interval = calling_intervals.key[0]
    (any_bad_intervals, chrs_represented) = calling_intervals.aggregate(
        (hl.agg.any(interval.start.contig != interval.end.contig), hl.agg.collect_as_set(interval.start.contig)))
    if any_bad_intervals:
        raise ValueError(
            "'impute_sex_chromosome_ploidy' does not support calling intervals that span chromosome boundaries")

    if len(rg.x_contigs) != 1:
        raise NotImplementedError(
            f"reference genome {rg.name!r} has multiple X contigs, this is not supported in 'impute_sex_chromosome_ploidy'"
        )
    if len(rg.y_contigs) != 1:
        raise NotImplementedError(
            f"reference genome {rg.name!r} has multiple Y contigs, this is not supported in 'impute_sex_chromosome_ploidy'"
        )

    kept_contig_filter = hl.array(chrs_represented).map(lambda x: hl.parse_locus_interval(x, reference_genome=rg))
    vds = VariantDataset(hl.filter_intervals(vds.reference_data, kept_contig_filter),
                         hl.filter_intervals(vds.variant_data, kept_contig_filter))

    mt = vds.variant_data
    calling_intervals = calling_intervals.annotate(interval_dup=interval)
    mt = mt.annotate_rows(interval=calling_intervals[mt.locus].interval_dup)
    mt = mt.filter_rows(hl.is_defined(mt.interval))
    coverage = mt.select_entries(sum_dp=mt.DP, interval_size=hl.is_defined(mt.DP))

    # checkpoint the coverage mt prior to returning
    checkpoint_path = tmp_prefix / 'sample_qc' / 'sex_coverage_precomputed.ht'
    if can_reuse(checkpoint_path):
        coverage = hl.read_matrix_table(str(checkpoint_path))
    else:
        coverage = coverage.checkpoint(str(checkpoint_path), overwrite=True)

    return coverage


def impute_sex(
    vds: hl.vds.VariantDataset,
    ht: hl.Table,
    tmp_prefix: Path,
) -> hl.Table:
    """
    Impute sex based on coverage.
    """
    checkpoint_path = tmp_prefix / 'sample_qc' / 'sex.ht'
    if can_reuse(str(checkpoint_path), overwrite=True):
        sex_ht = hl.read_table(str(checkpoint_path))
        return ht.annotate(**sex_ht[ht.s])

    # Load calling intervals
    seq_type = get_config()['workflow']['sequencing_type']
    calling_intervals_path = reference_path(f'broad/{seq_type}_calling_interval_lists')
    calling_intervals_ht = hl.import_locus_intervals(
        str(calling_intervals_path), reference_genome=genome_build()
    )
    logging.info('Calling intervals table:')
    calling_intervals_ht.describe()

    # clunky import due to dataproc execution
    from hail.vds.variant_dataset import VariantDataset

    # Pre-filter here and setting `variants_filter_lcr` and `variants_filter_segdup`
    # below to `False` to avoid the function calling gnomAD's `resources` module:
    for name in ['lcr_intervals_ht', 'seg_dup_intervals_ht']:
        interval_table = hl.read_table(str(reference_path(f'gnomad/{name}')))
        if interval_table.count() > 0:
            # remove all rows where the locus falls within a defined interval
            tmp_variant_data = vds.variant_data.filter_rows(
                hl.is_defined(interval_table[vds.variant_data.locus]), keep=False
            )
            vds = VariantDataset(
                reference_data=vds.reference_data, variant_data=tmp_variant_data
            ).checkpoint(str(tmp_prefix / f'{name}_checkpoint.vds'))
            logging.info(f'count post {name} filter:{vds.variant_data.count()}')

    # generate sex coverage mt, checkpointed
    coverage_mt = generate_sex_coverage_mt(vds, calling_intervals_ht, tmp_prefix)

    # Infer sex (adds row fields: is_female, var_data_chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        tmp_prefix=str(tmp_prefix / 'annotate_sex'),
        overwrite=not get_config()['workflow'].get('check_intermediates'),
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
        variants_only_x_ploidy=True,
        variants_only_y_ploidy=False,
        variants_filter_lcr=False,  # already filtered above
        variants_filter_segdup=False,  # already filtered above
        variants_filter_decoy=False,
        coverage_mt=coverage_mt
    )
    logging.info('Sex table:')
    sex_ht.describe()
    sex_ht = sex_ht.transmute(
        impute_sex_stats=hl.struct(
            f_stat=sex_ht.f_stat,
            n_called=sex_ht.n_called,
            expected_homs=sex_ht.expected_homs,
            observed_homs=sex_ht.observed_homs,
        )
    )
    sex_ht = sex_ht.checkpoint(str(checkpoint_path), overwrite=True)
    return sex_ht


def add_soft_filters(ht: hl.Table) -> hl.Table:
    """
    Uses the sex imputation results, variant sample qc, and input QC metrics
    to populate "filters" field to the sample table.
    """
    ht = ht.annotate(filters=hl.empty_set(hl.tstr))

    # Helper function to add filters into the `hard_filters` set
    def add_filter(ht_, expr, name):
        return ht_.annotate(
            filters=hl.if_else(
                expr & hl.is_defined(expr), ht_.filters.add(name), ht_.filters
            )
        )

    # Remove samples with ambiguous sex assignments
    ht = add_filter(ht, ht.sex_karyotype == 'ambiguous', 'ambiguous_sex')
    ht = add_filter(
        ht,
        ~hl.set({'ambiguous', 'XX', 'XY'}).contains(ht.sex_karyotype),
        'sex_aneuploidy',
    )

    cutoffs = get_config()['large_cohort']['sample_qc_cutoffs']
    ht = ht.annotate_globals(hard_filter_cutoffs=hl.struct(**cutoffs))

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    ht = add_filter(
        ht,
        ht.var_data_chr20_mean_dp < cutoffs['min_coverage'],
        'low_coverage',
    )

    # Remove extreme raw bi-allelic sample QC outliers
    ht = add_filter(
        ht,
        (
            (ht.sample_qc.n_snp > cutoffs['max_n_snps'])
            | (ht.sample_qc.n_snp < cutoffs['min_n_snps'])
            | (ht.sample_qc.n_singleton > cutoffs['max_n_singletons'])
            | (ht.sample_qc.r_het_hom_var > cutoffs['max_r_het_hom'])
        ),
        'bad_sample_qc_metrics',
    )
    ht = ht.annotate(filtered=hl.len(ht.filters) > 0)
    logging.info('Table with filters:')
    ht.describe()
    return ht
