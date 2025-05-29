import os
from lib2to3.pgen2 import driver

import click

import hail as hl

from cpg_utils import Path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch, init_batch, reference_path
from gnomad.sample_qc.pipeline import get_qc_mt

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
    '--capture-region-bed-files',
    help='Names of capture region bed files to use for exomes. Names must be keys in references dictionary.',
    type=str,
    multiple=True,
)
@click.command()
def main(vds_path: str, exomes: bool, bed_files: list[str]):
    print(f'Input vds_path: {vds_path}')
    pruned_variant_table_path = output_path('pruned_variants_exome.ht', 'default')
    print('Will be writing to pruned_variant_table_path:', pruned_variant_table_path)
    # Initialise batch
    init_batch(
        worker_memory='highmem',
        driver_memory='highmem',
        driver_cores=4,
    )
    # exomes
    if exomes:
        if not bed_files:
            raise ValueError('If --exomes is set, you must provide at least one --capture-region-bed-files')

        # Read in capture region bed files
        intervals = []
        for bed_file in bed_files:
            bed_path = reference_path(f'{bed_file}')
            if not os.path.exists(bed_path):
                raise FileNotFoundError(f'Capture region bed file {bed_path} does not exist')
            intervals.append(hl.import_bed(str(bed_path), reference_genome='GRCh38'))

        # Intersect the intervals
        intervals = hl.intersection(*intervals)
    # if exomes read in capture region bed files
    # generate hl.Interval objects for each capture region in bed files
    # pass list of hl.Interval objects to hl.vds.read_vds(vds_path, intervals=intervals)
    # the intervals need to be the INTERSECTION not the UNION of the capture regions
    vds = hl.vds.read_vds(str(vds_path))

    print('Splitting multi-allelic sites')
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    print('Done splitting multi-allelic sites')

    print('Densifying VDS')
    hgdp_1kg = hl.vds.to_dense_mt(vds)
    print('Done densifying VDS')

    # Run variant QC
    print('Running variant QC')
    # choose variants based off of gnomAD v3 parameters
    # Inbreeding coefficient > -0.25 (no excess of heterozygotes)
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
        & (hgdp_1kg.IB.f_stat > -0.25),
    )
    print('Done filtering using gnomAD v3 parameters')

    print('Filtering to exome intervals')
    exome_intervals = hl.import_locus_intervals(str(reference_path('broad/exome_calling_interval_lists')))
    hgdp_1kg_exome = hgdp_1kg.filter_rows(hl.is_defined(exome_intervals[hgdp_1kg.locus]))
    print('Done filtering to exome intervals')

    # downsize input variants for ld_prune
    # otherwise, persisting the pruned_variant_table will cause
    # script to fail. See https://github.com/populationgenomics/ancestry/pull/79
    print('Writing sites table pre-LD pruning')
    checkpoint_path = output_path('hgdp1_1kg_exome_pre_pruning.mt', 'default')
    hgdp_1kg_exome = hgdp_1kg_exome.checkpoint(checkpoint_path, overwrite=True)
    print('Done writing sites table pre-LD pruning')

    # nrows = hgdp_1kg_exome.count_rows()
    # print(f'hgdp_1kg_exome.count_rows() = {nrows}')
    # hgdp_1kg_exome = hgdp_1kg_exome.sample_rows(
    #     NUM_ROWS_BEFORE_LD_PRUNE / nrows,
    #     seed=12345,
    # )

    # as per gnomAD, LD-prune variants with a cutoff of r2 = 0.1
    print('Pruning sites table')
    pruned_variant_table = hl.ld_prune(
        hgdp_1kg_exome.GT,
        r2=0.1,
        bp_window_size=500000,
    )
    print('Done pruning sites table')

    print('Repartitioning sites table')
    # repartition table after pruning
    pruned_variant_table = pruned_variant_table.repartition(100, shuffle=False)
    print('Number of variants in pruned_variant_table:', pruned_variant_table.count())
    print('Done repartitioning sites table')
    # pruned_variant_table_path = output_path('pruned_variants_exome.ht', 'tmp')
    print(f'Writing sites table to {pruned_variant_table_path}')
    pruned_variant_table.write(pruned_variant_table_path, overwrite=True)
    print('Done writing sites table')


if __name__ == '__main__':
    main()
