#!/usr/bin/env python3

"""
Takes a MatrixTable and two output paths
Writes two representations - region filtered MatrixTable, and a sites-only VCF representation
VCF is exported per-interval, with a separate header file, to be concatenated later

This script also takes a BED file as input; output contains only the variants that overlap with the BED file

All existing info fields are dropped, and replaced with just the callset AC / AN / AF

This removes any Filtered variants, due to an issue between the header and content
 - when we apply VQSR annotations we pull in Filters, but we don't pull the corresponding header lines
 - this means that the MT contains rows which aren't explained in the header, causing some tools to fail
 - for now it's easier to just remove these rows - reconsider if we use this properly
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import config_retrieve, genome_build
from cpg_utils.hail_batch import init_batch
from cpg_workflows.utils import get_logger

# VQSR filters are added to the VQSR VCF, but we don't retain these header lines when we copy the VQSR info over
# to the variant MatrixTable during AnnotateCohort.
# If we export VCFs here without these header lines, there's a possibility that downstream tools won't be able to
# process them correctly, as they won't be able to find the VQSR filters in the VCF Header.
# We can either pair this workflow more tightly to the combiner workflow, and pass in both a MatrixTable and a VQSR VCF
# Or we can just bundle a dictionary here of all the VQSR headers we could be applying.
VQSR_FILTERS = {
    'filter': {
        "VQSRTrancheINDEL99.00to99.50": {
            "Description": "Truth sensitivity tranche level for INDEL model at VQS Lod: -1.4652 <= x < -0.6489",
        },
        "VQSRTrancheINDEL99.50to99.90": {
            "Description": "Truth sensitivity tranche level for INDEL model at VQS Lod: -8.3914 <= x < -1.4652",
        },
        "VQSRTrancheINDEL99.90to99.95": {
            "Description": "Truth sensitivity tranche level for INDEL model at VQS Lod: -20.9224 <= x < -8.3914",
        },
        "VQSRTrancheINDEL99.95to100.00": {
            "Description": "Truth sensitivity tranche level for INDEL model at VQS Lod: -39995.8675 <= x < -20.9224",
        },
        "VQSRTrancheINDEL99.95to100.00+": {
            "Description": "Truth sensitivity tranche level for INDEL model at VQS Lod < -39995.8675",
        },
        "VQSRTrancheSNP99.00to99.90+": {
            "Description": "Truth sensitivity tranche level for SNP model at VQS Lod < -10.0",
        },
        "VQSRTrancheSNP99.90to100.00": {
            "Description": "Truth sensitivity tranche level for SNP model at VQS Lod: -10.0 <= x < -4.37",
        },
        "VQSRTrancheSNP99.90to100.00+": {
            "Description": "Truth sensitivity tranche level for SNP model at VQS Lod < -10.0",
        },
    },
}


def main(
    mt_path: str,
    output_mt: str,
    output_sites_only: str,
    bed: str | None,
) -> None:
    """

    Args:
        mt_path (str):
        output_mt (str): write region-filtered MatrixTable, stripped of INFO fields
        output_sites_only (str): write a per-partition sites-only VCF directory to this location
        bed (str): Region BED file
    """

    init_batch(
        worker_memory=config_retrieve(['combiner', 'worker_memory'], 'highmem'),
        driver_memory=config_retrieve(['combiner', 'driver_memory'], 'highmem'),
        driver_cores=config_retrieve(['combiner', 'driver_cores'], 2),
    )

    # read the dense MT and obtain the sites-only HT
    mt = hl.read_matrix_table(mt_path)

    if bed:
        # remote-read of the BED file, skipping any contigs not in the reference genome
        # the Ensembl data wasn't filtered to remove non-standard contigs
        limited_region = hl.import_bed(bed, reference_genome=genome_build(), skip_invalid_intervals=True)
        # filter to overlaps with the BED file
        mt = mt.filter_rows(hl.is_defined(limited_region[mt.locus]))

    # replace the existing INFO block to just have AC/AN/AF - no other carry-over
    # this is based on the structure we already achieved in annotate_cohort
    # drop every other field (scrapping all the previously generated VEP/ClinVar/etc)
    # these fields are mangled a little in AnnotateCohort - revisit that
    mt = mt.select_rows(
        info=hl.struct(
            AF=mt.info.AF,
            AN=mt.info.AN,
            AC=mt.info.AC,
        ),
        rsid=mt.rsid,
        filters=mt.filters,
    )

    mt.describe()

    mt.write(
        output_mt,
        overwrite=True,
    )

    # now read that location for speed, and write the sites-only VCF
    # keep partitions consistent
    sites_only_ht = hl.read_matrix_table(output_mt).rows()

    get_logger().info('Writing sites-only VCF in fragments, header-per-shard')
    hl.export_vcf(sites_only_ht, output_sites_only, tabix=True, parallel='separate_header', metadata=VQSR_FILTERS)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to the input MT',
        required=True,
    )
    parser.add_argument(
        '--output_mt',
        help='Path to write the resulting MatrixTable',
        required=True,
    )
    parser.add_argument(
        '--output_sites_only',
        help='Specify an output path for a sites-only VCF, or None',
    )
    parser.add_argument(
        '--bed',
        help='Region BED file',
    )
    args = parser.parse_args()
    main(mt_path=args.input, output_mt=args.output_mt, output_sites_only=args.output_sites_only, bed=args.bed)


if __name__ == '__main__':
    cli_main()
