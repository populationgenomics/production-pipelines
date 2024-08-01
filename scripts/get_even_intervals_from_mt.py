"""
Read a matrixtable and split it into even intervals of rows.
Use the contigs at the start of each interval to define the
intervals saved in the output interval_list file.

Hard breaks at centromeres.
"""

import argparse
import json

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import dataset_path, reference_path
from cpg_utils.hail_batch import init_batch

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def get_intervals_for_chr(
    mt: hl.MatrixTable,
    interval_size: int,
    telomere_positions: dict[str, tuple[int, int]],
    centromere_position: tuple[int, int],
) -> list[hl.MatrixTable]:
    """
    Evenly split the matrixtable into intervals for variants in a single chromosome.
    Hard breaks at telomeres and centromeres.
    """

    # The first valid region for genotyping begins after the first telomere and ends before the centromere
    genotyping_region_1_start = telomere_positions['first'][1] + 1
    genotyping_region_1_end = centromere_position[0] - 1
    mt_region_1 = mt.filter_rows(
        (mt.locus.position >= genotyping_region_1_start) & (mt.locus.position < genotyping_region_1_end),
    )
    mt_region_1 = mt_region_1.add_row_index(name='row_idx')
    mt_region_1 = mt_region_1.rows().select('row_idx', 'locus')

    # The second valid region for genotyping starts after the centromere and ends before the second telomere
    genotyping_region_2_start = centromere_position[1] + 1
    genotyping_region_2_end = telomere_positions['second'][0] - 1
    mt_region_2 = mt.filter_rows(
        (mt.locus.position >= genotyping_region_2_start) & (mt.locus.position < genotyping_region_2_end),
    )
    mt_region_2 = mt_region_2.add_row_index(name='row_idx')
    mt_region_2 = mt_region_2.rows().select('row_idx', 'locus')

    # Iterate through the first genotyping region mt in intervals of interval_size, using the row_idx
    intervals = []
    for mt_region in [mt_region_1, mt_region_2]:
        for i in range(0, mt_region.count_rows(), interval_size):
            # Get the interval
            interval = mt_region.filter_rows((mt_region.row_idx >= i) & (mt_region.row_idx < i + interval_size))

            # Add the interval to the intervals list
            intervals.append(interval)

    return intervals


def get_even_intervals_from_mt(
    mt: hl.MatrixTable,
    n_intervals: int,
    tc_intervals_filepath: str | None,
    output_intervals_path: str,
):
    """
    Get even intervals from a matrixtable.
    """
    # Get the number of rows in the matrixtable
    n_rows = mt.count_rows()

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals

    # Get the telomere and centromere positions
    telomere_positions, centromere_positions = get_telomere_and_centromere_start_end_positions(tc_intervals_filepath)

    intervals = {}
    for chrom in CHROMS:
        chrom_intervals = get_intervals_for_chr(
            mt=mt.filter_rows(mt.locus.contig == chrom),
            interval_size=interval_size,
            telomere_positions=telomere_positions[chrom],
            centromere_position=centromere_positions[chrom],
        )
        intervals[chrom] = chrom_intervals

    # Evaluate the intervals
    interval_positions_by_chrom = write_intervals_json(intervals, output_intervals_path)

    return interval_positions_by_chrom


def write_intervals_json(intervals: dict[str, list[hl.MatrixTable]], output_path: str):
    """
    Evaluate the list of MatrixTable intervals and save them to a json.
    """
    interval_positions: dict[str, list[tuple[int, int]]] = {}
    for chrom, chrom_intervals in intervals.items():
        interval_positions[chrom] = []
        for i, interval in enumerate(chrom_intervals):
            print(f'Chrom {chrom}, interval {i}: {interval.count_rows()} rows')
            # Get the start and end position of the interval
            positions = interval.locus.position.collect()
            start_pos = positions[0]
            end_pos = positions[-1]
            print(f'Interval start: {start_pos}, end: {end_pos}')
            interval_positions[chrom].append((start_pos, end_pos))

    # Save the interval positions to a file
    with open(output_path, 'w') as f:
        json.dump(interval_positions, f)

    print(f'Intervals saved to {output_path}')
    return interval_positions


def get_telomere_and_centromere_start_end_positions(
    tc_intervals_filepath: str | None,
) -> tuple[dict[str, dict[str, tuple[int, int]]], dict[str, tuple[int, int]]]:
    """
    Get the start and end positions of the telomeres and centromeres in the human genome.
    """
    tc_intervals_path = (
        to_path(tc_intervals_filepath)
        if tc_intervals_filepath
        else to_path(
            reference_path('hg38_telomeres_and_centromeres_intervals/interval_list'),
        )
    )

    with open(tc_intervals_path) as f:
        # Read the intervals from the source file, skip lines starting with '@'
        tc_intervals = [line.strip().split() for line in f if not line.startswith('@')]

    # Collect the start and end positions of the centromeres and telomeres by chromosome
    centromeres_by_chrom = {}
    telomeres_by_chrom: dict[str, dict[str, tuple[int, int]]] = {}
    for chrom, start, end, _, region in tc_intervals:
        if chrom not in telomeres_by_chrom:
            telomeres_by_chrom[chrom] = {}
        if 'centromere' in region:
            centromeres_by_chrom[chrom] = (int(start), int(end))
        if 'telomere' in region:
            if start == '1':
                # If the start position is '1', the telomere is at the start of the chromosome
                telomeres_by_chrom[chrom].update({'first': (int(start), int(end))})
            else:
                # If the start position is not '1', the telomere is at the end of the chromosome
                telomeres_by_chrom[chrom].update({'second': (int(start), int(end))})

    return telomeres_by_chrom, centromeres_by_chrom


def main(args):
    """
    Run the script.
    """
    output_intervals_path = dataset_path('derived_intervals/hg38_even_intervals.json')
    init_batch()
    mt = hl.read_matrix_table(args.input_mt)
    get_even_intervals_from_mt(mt, args.n_intervals, args.tc_intervals_filepath, output_intervals_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-mt', help='Input matrixtable', required=True)
    parser.add_argument(
        '--tc-intervals-filepath',
        help='Telomere and centromere interval_list/bed file path',
        required=False,
        default=None,
    )
    parser.add_argument(
        '--n-intervals',
        help='Number of intervals to split the matrixtable into',
        type=int,
        required=True,
    )
    args = parser.parse_args()
    main(args)
