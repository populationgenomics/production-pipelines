"""
Read a matrixtable and split it into even intervals of rows.
Use the contigs at the start of each interval to define the
intervals saved in the output interval_list file.

Hard breaks at centromeres.
"""

import argparse
import json
import logging

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import dataset_path, reference_path
from cpg_utils.hail_batch import init_batch

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def get_telomere_and_centromere_start_end_positions(
    tc_intervals_filepath: str | None,
) -> tuple[dict[str, dict[str, tuple[int, int]]], dict[str, tuple[int, int]]]:
    """
    Get the start and end positions of the telomeres and centromeres in the human genome.
    Each chromosome has two telomeres and one centromere.
    Returns:
        telomeres_by_chrom:   {'chr1': {'first': (1, 10000), 'second': (248946423, 248956422)}, ...}
        centromeres_by_chrom: {'chr1': (122026460, 124932724), ...}
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


def evaluate_ht_positions(ht: hl.Table) -> tuple[list[str], list[int], list[int]]:
    """
    Call hl.Table.aggregate() to collect the positions of the intervals.
    Returns:
        chroms:     ['chr1', 'chr1', ...]
        positions:  [start_1, stop_1, start_2, stop_2, ...]
        idx:        [X1, X2, X3, X4, ...]
    """
    struct = hl.struct(
        chrom=ht.locus.contig,
        pos=ht.locus.position,
        idx=ht.global_row_idx,
    )
    result = ht.aggregate(hl.agg.collect(struct))
    chroms = [r.chrom for r in result]
    positions = [r.pos for r in result]
    idx = [r.idx for r in result]

    return chroms, positions, idx


def get_intervals_from_intervals_ht(
    intervals_ht: hl.Table,
) -> dict[str, list[tuple[int, int]]]:
    """
    Evaluate the start and end of each interval and return a dictionary of the intervals by chromosome.
    """
    chroms, positions, idx = evaluate_ht_positions(intervals_ht)

    # Create dictionaries of the interval start and stop positions by chromosome
    interval_starts_by_chrom = {}
    interval_ends_by_chrom = {}
    for counter, (chr, pos, i) in enumerate(zip(chroms, positions, idx)):
        if chr not in interval_starts_by_chrom:
            interval_starts_by_chrom[chr] = []
            interval_ends_by_chrom[chr] = []
        # odd iterations are the starts of the intervals, evens are the ends
        if counter % 2 != 0:
            interval_starts_by_chrom[chr].append(pos)
            logging.info(f'{chr} :: start :: {pos} :: idx :: {i}')
        else:
            interval_ends_by_chrom[chr].append(pos)
            logging.info(f'{chr} :: stop  :: {pos} :: idx :: {i}')

    intervals_by_chrom = {}
    for chrom in CHROMS:
        intervals_by_chrom[chrom] = []
        for start, end in zip(interval_starts_by_chrom[chrom], interval_ends_by_chrom[chrom]):
            intervals_by_chrom[chrom].append((start, end))

    return intervals_by_chrom


def get_intervals_from_ht(
    ht: hl.Table,
    n_intervals: int,
    tc_intervals_filepath: str | None,
    output_intervals_path: str,
):
    """
    Read the matrixtable rows (ht), and split them into n_intervals of 
    approximately equal sized rows.

    Subsets the ht by locus.contig to get the intervals for each chromosome. 
    Hard exclude the telomeric and centromeric regions of each chromosome by 
    reading the telomere and centromere positions from the tc_intervals_filepath.

    Return the start and end positions of the intervals in a json file.
    """
    # Get the number of rows in the table
    n_rows = ht.count()
    logging.info(f'Number of rows in the table: {n_rows}')

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals
    logging.info(f'Interval size: {interval_size}')

    # Get the telomere and centromere positions
    telomere_positions, centromere_positions = get_telomere_and_centromere_start_end_positions(tc_intervals_filepath)

    # Get the intervals for each chromosome 
    intervals: list[hl.Table] = []
    for chrom in CHROMS:
        intervals.extend(
            get_ht_intervals_for_chrom(
                ht=ht.filter(ht.locus.contig == chrom),
                interval_size=interval_size,
                telomere_positions=telomere_positions[chrom],
                centromere_position=centromere_positions[chrom],
            ),
        )

    # Stack the hail Table intervals into a single Table
    intervals_ht = hl.Table.union(*intervals)
    intervals_by_chrom = get_intervals_from_intervals_ht(intervals_ht, n_intervals, output_intervals_path)

    with to_path(output_intervals_path).open('w') as f:
        json.dump(intervals_by_chrom, f)
    logging.info(f'Intervals saved to {output_intervals_path}')

    return intervals_by_chrom


def filter_telomeres_and_centromeres(ht, telomere_positions, centromere_position):
    """
    Filter the table to exclude the telomeric and centromeric regions of the chromosome.
    """
    # Filter the table to exclude the telomeric regions
    ht = ht.filter(
        (ht.locus.position < telomere_positions['first'][0]) | (ht.locus.position > telomere_positions['first'][1]),
    )
    ht = ht.filter(
        (ht.locus.position < telomere_positions['second'][0]) | (ht.locus.position > telomere_positions['second'][1]),
    )

    # Filter the table to exclude the centromeric region
    ht = ht.filter(
        (ht.locus.position < centromere_position[0]) | (ht.locus.position > centromere_position[1]),
    )

    return ht


def get_ht_intervals_for_chrom(
    ht: hl.Table,
    interval_size: int,
    telomere_positions: dict[str, tuple[int, int]],
    centromere_position: tuple[int, int],
) -> list[hl.Table]:
    """
    Evenly split the chromosome subset of the table into intervals of interval_size.
    Hard exclude the telomeric and centromeric regions of the chromosome.

    Return a list of hail Tables, each containing the start or end row of an interval.
    """
    ht = filter_telomeres_and_centromeres(ht, telomere_positions, centromere_position)

    ht_region_1 = ht.filter(ht.locus.position < centromere_position[0]).add_index(name='row_idx')
    ht_region_2 = ht.filter(ht.locus.position > centromere_position[1]).add_index(name='row_idx')

    # Iterate through both genotyping regions in intervals of interval_size, using the row_idx
    intervals = []
    for ht_region in [ht_region_1, ht_region_2]:
        for i in range(0, ht_region.count(), interval_size):
            # Get the interval start and end rows
            interval_start = ht_region.filter(ht_region.row_idx == i)
            # The last interval may be smaller than interval_size
            interval_end = ht_region.filter(ht_region.row_idx == min((i + interval_size - 1), ht_region.count() - 1))

            # Add the interval start and end rows to the intervals list
            intervals.extend([interval_start, interval_end])

    return intervals


def main(args):
    """
    Run the script.
    """
    init_batch(driver_memory='16G')
    mt = hl.read_matrix_table(args.input_mt)
    logging.info('Read matrixtable')
    ht = mt.rows()
    ht = ht.select().select_globals()
    ht = ht.add_index(name='global_row_idx')
    ht = ht.key_by('locus', 'alleles', 'global_row_idx')
    logging.info('Re-keyed table')
    intervals = get_intervals_from_ht(ht, args.n_intervals, args.tc_intervals_filepath, args.output_intervals_path)
    logging.info('Calculated intervals')


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
    parser.add_argument(
        '--output-intervals-path',
        help='Path to save the intervals file',
        required=False,
        default=dataset_path('derived_intervals/hg38_even_intervals.json'),
    )
    args = parser.parse_args()
    main(args)
