"""
Impute sex.
Add soft filters for samples.
"""

import logging

import hail as hl
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, genome_build
from cpg_workflows.inputs import get_cohort
from cpg_workflows.utils import can_reuse
from gnomad.sample_qc.pipeline import annotate_sex


def run(
    vds_path: Path,
    out_sample_qc_ht_path: Path,
    tmp_prefix: Path,
):
    if can_reuse(out_sample_qc_ht_path):
        return hl.read_table(str(out_sample_qc_ht_path))

    ht = initialise_sample_table()

    vds = hl.vds.read_vds(str(vds_path))

    # Remove centromeres and telomeres:
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    if tel_cent_ht.count() > 0:
        vds = hl.vds.filter_intervals(vds, tel_cent_ht, keep=False)

    # Run Hail sample-QC stats:
    sqc_ht_path = tmp_prefix / 'sample_qc.ht'
    if can_reuse(sqc_ht_path):
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
    pop_meta_field = get_config()['large_cohort'].get('pop_meta_field')
    entries = [
        {
            's': s.id,
            'external_id': s.external_id,
            'dataset': s.dataset.name,
            'gvcf': str(s.gvcf.path),
            'sex': s.pedigree.sex.value,
            'pop': s.meta.get(pop_meta_field) if pop_meta_field else None,
        }
        for s in get_cohort().get_samples()
        if s.gvcf
    ]
    if not entries:
        raise ValueError('No samples with GVCFs found')
    t = 'array<struct{s: str, external_id: str, dataset: str, gvcf: str, sex: int, pop: str}>'
    ht = hl.Table.parallelize(hl.literal(entries, t), key='s')
    return ht


def impute_sex(
    vds: hl.vds.VariantDataset,
    ht: hl.Table,
    tmp_prefix: Path,
) -> hl.Table:
    """
    Impute sex based on coverage.
    """
    checkpoint_path = tmp_prefix / 'sample_qc' / 'sex.ht'
    if can_reuse(str(checkpoint_path)):
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

    # Pre-filter here and setting `variants_filter_lcr` and `variants_filter_segdup`
    # below to `False` to avoid the function calling gnomAD's `resources` module:
    for name in ['lcr_intervals_ht', 'seg_dup_intervals_ht']:
        ht = hl.read_table(str(reference_path(f'gnomad/{name}')))
        if ht.count() > 0:
            vds = hl.vds.filter_intervals(vds, ht, keep=False)

    # Infer sex (adds row fields: is_female, var_data_chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        tmp_prefix=str(tmp_prefix / 'annotate_sex'),
        overwrite=not get_config()['workflow'].get('check_intermediates'),
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
        variants_only_x_ploidy=True,
        variants_only_y_ploidy=True,
        variants_filter_lcr=False,  # already filtered above
        variants_filter_segdup=False,  # already filtered above
        variants_filter_decoy=False,
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
