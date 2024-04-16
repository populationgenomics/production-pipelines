import click

import hail as hl

from cpg_utils import Path
from cpg_utils.hail_batch import get_batch, init_batch, output_path, reference_path
from gnomad.sample_qc.pipeline import get_qc_mt


def sites_table(vds_path: str, gcs_output_path: str) -> hl.MatrixTable:
    """
    Select variants within high-quality intervals according to gnomAD v4 criteria.

    Input:
        - Combined VDS with all datasets (exomes and genomes)
            - Exomes should be filtered to high-quality intervals prior to running this script
    Output:
        - sites_table: MatrixTable containing high-quality sites
    """
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
    print('Variant count before split:', vds.variant_data.count())
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    print('Variant count after split:', vds.variant_data.count())
    print('Done splitting multi-allelic sites')

    print('Densifying VDS')
    print('Variant count before densify:', vds.variant_data.count())
    mt = hl.vds.to_dense_mt(vds)
    print('Variant count after densify:', mt.count())
    mt.show(n_cols=10)
    print('Done densifying VDS')

    print('Generating sites table')
    sites_table = get_qc_mt(
        mt,
        min_af=0.0001,
        min_callrate=0.95,
        filter_lcr=False,  # Already filtered above
        filter_segdup=False,  # Already filtered above
        checkpoint_path=gcs_output_path,
    )
    print('Done generating sites table')

    print(f'Writing sites table to {gcs_output_path}')
    sites_table.checkpoint(gcs_output_path, overwrite=True)
    print('Done writing sites table')


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

    gcs_output_path = 'gs://cpg-thousand-genomes-test/exome/sites_table/'

    print('Generating sites table')
    sites_table = get_qc_mt(
        mt,
        min_af=0.0001,
        min_callrate=0.95,
        filter_lcr=False,  # Already filtered above
        filter_segdup=False,  # Already filtered above
        filter_decoy=False,  # gnomad.resources.grch38.reference_data has no attribute 'decoy_intervals'
        checkpoint_path=gcs_output_path,
    )
    print('Done generating sites table')

    # Initialise job

    print('here1 and vds_path:', vds_path, type(vds_path))
    print('here2')
    print('here3 and vds_path:', vds_path, type(vds_path))
    print('here3 and output_path:', gcs_output_path, type(gcs_output_path))

    print(f'Writing sites table to {gcs_output_path}')
    sites_table.checkpoint(gcs_output_path, overwrite=True)
    print('Done writing sites table')


if __name__ == '__main__':
    main()
