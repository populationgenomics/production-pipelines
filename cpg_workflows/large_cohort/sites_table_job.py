#!/usr/bin/env python3

import click

import hail as hl

from cpg_utils import Path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import genome_build, init_batch

NUM_ROWS_BEFORE_LD_PRUNE = 200000


@click.option(
    '--vds-path',
    help='Path to the VDS file',
    type=str,
)
@click.option(
    '--exomes',
    help='Whether to process exomes. If true, will read in capture region bed files.',
    is_flag=True,
    default=False,
)
@click.option(
    '--intersected-bed-file',
    help='A bed file that contains the intersection of all capture regions. '
    'If --exomes is set, this must be provided.',
    type=str,
)
@click.option(
    '--vqsr-table-path',
    help='Path to the VQSR table to use for filtering variants.',
    type=str,
)
@click.option(
    '--sites-table-outpath',
    help='Path to write the output sites table.',
    type=str,
)
@click.command()
def main(
    vds_path: str,
    exomes: bool,
    vqsr_table_path: str,
    sites_table_outpath: str,
    intersected_bed_file: list[str] | None = None,
):
    print(f'Input vds_path: {vds_path}')

    print('Will be writing to pruned_variant_table_path:', sites_table_outpath)

    # Initialise batch
    init_batch(
        worker_memory='highmem',
        driver_memory='highmem',
        driver_cores=4,
    )
    vqsr_table = hl.read_table(vqsr_table_path)

    # exomes
    if exomes:
        if not intersected_bed_file:
            raise ValueError('If --exomes is set, you must provide at least one --capture-region-bed-files')

        # Read in capture region bed files
        capture_interval_ht: hl.Table = hl.import_bed(str(intersected_bed_file), reference_genome=genome_build())

        # Generate list of intervals
        intervals: list[hl.Interval] = capture_interval_ht.interval.collect()

    vds = hl.vds.read_vds(str(vds_path), intervals=intervals if exomes else None)

    # Filter to variant sites that pass VQSR
    passed_variants = vqsr_table.filter(vqsr_table.info.AS_FilterStatus == 'PASS')
    vds = hl.vds.filter_variants(vds, passed_variants)

    print('Splitting multi-allelic sites')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    print('Done splitting multi-allelic sites')

    print('Densifying VDS')
    hgdp_1kg = hl.vds.to_dense_mt(vds)
    print('Done densifying VDS')

    # Run variant QC
    print('Running variant QC')
    # choose variants based off of gnomAD v3 parameters
    # Inbreeding coefficient > -0.80 (no excess of heterozygotes)
    # Must be single nucleotide variants that are autosomal (i.e., no sex), and bi-allelic
    # Have an allele frequency above 1% (note deviation from gnomAD, which is 0.1%)
    # Have a call rate above 99%
    hgdp_1kg = hl.variant_qc(hgdp_1kg)
    print('Done running variant QC')

    print('Generating sites table')
    print('Filtering using gnomAD v3 parameters')
    hgdp_1kg = hgdp_1kg.annotate_rows(
        IB=hl.agg.inbreeding(
            hgdp_1kg.GT,
            hgdp_1kg.variant_qc.AF[1],
        ),
    )
    hgdp_1kg = hgdp_1kg.filter_rows(
        (hl.len(hgdp_1kg.alleles) == 2)
        & (hgdp_1kg.locus.in_autosome())
        & (hgdp_1kg.variant_qc.AF[1] > 0.01)
        & (hgdp_1kg.variant_qc.call_rate > 0.99)
        & (hgdp_1kg.IB.f_stat > -0.80),
    )
    print('Done filtering using gnomAD v3 parameters')

    # downsize input variants for ld_prune
    # otherwise, persisting the pruned_variant_table will cause
    # script to fail. See https://github.com/populationgenomics/ancestry/pull/79
    print('Writing sites table pre-LD pruning')
    checkpoint_path = output_path('hgdp_1kg_exome_pre_pruning.mt', 'default')
    hgdp_1kg = hgdp_1kg.checkpoint(checkpoint_path, overwrite=True)
    print('Done writing sites table pre-LD pruning')

    # nrows = hgdp_1kg.count_rows()
    # print(f'hgdp_1kg.count_rows() = {nrows}')
    # hgdp_1kg = hgdp_1kg.sample_rows(
    #     NUM_ROWS_BEFORE_LD_PRUNE / nrows,
    #     seed=12345,
    # )

    # as per gnomAD, LD-prune variants with a cutoff of r2 = 0.1
    print('Pruning sites table')
    pruned_variant_table = hl.ld_prune(
        hgdp_1kg.GT,
        r2=0.1,
        bp_window_size=500000,
    )
    print('Done pruning sites table')

    print('Repartitioning sites table')
    # repartition table after pruning
    pruned_variant_table = pruned_variant_table.repartition(100, shuffle=False)
    print('Number of variants in pruned_variant_table:', pruned_variant_table.count())
    print('Done repartitioning sites table')

    print(f'Writing sites table to {sites_table_outpath}')
    pruned_variant_table.write(sites_table_outpath, overwrite=True)
    print('Done writing sites table')


if __name__ == '__main__':
    main()
