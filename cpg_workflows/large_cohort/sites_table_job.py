import os

import click

import hail as hl

from cpg_utils import Path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch, init_batch, reference_path
from gnomad.sample_qc.pipeline import get_qc_mt


@click.option(
    '--vds-path',
    help='Path to the VDS file',
    type=str,
)
@click.command()
def main(vds_path):
    # Initialise batch
    init_batch()

    vds = hl.vds.read_vds(str(vds_path))

    # Pre-filtering
    print('Filtering centromeres and telomeres')
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    if tel_cent_ht.count() > 0:
        vds = hl.vds.filter_intervals(vds, tel_cent_ht, keep=False)
    print('Done filtering centromeres and telomeres')

    print('Filtering to autosomes')
    vds = hl.vds.filter_chromosomes(vds, keep=[f'chr{chrom}' for chrom in range(1, 23)])
    print('Done filtering to autosomes')

    print('Filtering low complexity regions and segmental duplications')
    for name in ['lcr_intervals_ht', 'seg_dup_intervals_ht']:
        interval_table = hl.read_table(str(reference_path(f'gnomad/{name}')))
        if interval_table.count() > 0:
            # remove all rows where the locus falls within a defined interval
            tmp_variant_data = vds.variant_data.filter_rows(
                hl.is_defined(interval_table[vds.variant_data.locus]),
                keep=False,
            )
            vds = hl.vds.VariantDataset(reference_data=vds.reference_data, variant_data=tmp_variant_data)
    print('Done filtering low complexity regions and segmental duplications')

    print('Splitting multi-allelic sites')
    # print('Variant count before split:', vds.variant_data.count())
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    # print('Variant count after split:', vds.variant_data.count())
    print('Done splitting multi-allelic sites')

    print('Densifying VDS')
    # print('Variant count before densify:', vds.variant_data.count())
    mt = hl.vds.to_dense_mt(vds)
    # print('Variant count after densify:', mt.count())
    # mt.show(n_cols=10)
    print('Done densifying VDS')

    # Run variant QC
    print('Running variant QC')
    # Choosing variants based off of gnomAD v4 parameters
    # Inbreeding coefficient > -0.8
    # Bi-allelic, autosomal variants only
    # Allele frequency above 0.01%
    # Call rate above 95%
    # HWE p-value above 1e-8
    mt = hl.variant_qc(mt)
    print('Done running variant QC')

    print('Generating sites table')
    print('Filtering using gnomAD v4 parameters')
    mt = mt.annotate_rows(
        IB=hl.agg.inbreeding(
            mt.GT,
            mt.variant_qc.AF[1],
        ),
    )
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & (mt.locus.in_autosome())
        & (mt.variant_qc.AF[1] > 0.0001)
        & (mt.variant_qc.call_rate > 0.95)
        & (mt.IB.f_stat > -0.8)
        & (mt.variant_qc.p_value_hwe > 1e-8),
    )
    print('Done filtering using gnomAD v4 parameters')

    print('Writing sites table pre-LD pruning')
    checkpoint_path = output_path('sites_table.mt', 'tmp')
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    print('Done writing sites table pre-LD pruning')

    print('Pruning sites table')
    pruned_variant_table = hl.ld_prune(
        mt.GT,
        r2=0.1,
        bp_window_size=500000,
    )
    print('Done pruning sites table')

    print('Repartitioning sites table')
    # repartition table after pruning
    pruned_variant_table = pruned_variant_table.repartition(100, shuffle=False)
    print('Done repartitioning sites table')
    pruned_variant_table.checkpoint(output_path('final_sites_table.mt', 'default'), overwrite=True)
    print('Done writing sites table')


if __name__ == '__main__':
    main()
