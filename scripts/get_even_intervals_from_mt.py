"""
Read a matrixtable and split it into even intervals of rows.
Use the contigs at the start of each interval to define the
intervals saved in the output interval_list file.

Hard breaks at centromeres.
"""

import argparse
import json
import logging
import time
from typing import Generator, List, NamedTuple, Optional

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import dataset_path, reference_path
from cpg_utils.hail_batch import init_batch


class IntervalInfo(NamedTuple):
    start_idx: int
    end_idx: int
    start_pos: int
    end_pos: int
    start_contig: Optional[str]
    end_contig: Optional[str]

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


def get_even_intervals_from_mt_rows(
    ht: hl.Table,
    n_intervals: int,
    tc_intervals_filepath: str | None,
    output_intervals_path: str,
):
    """
    Read the matrixtable rows (ht) and split it into even intervals of rows by subsetting the rows by chromosome.
    Return the start and end positions of the intervals in a json file.
    """
    # Get the number of rows in the table
    n_rows = ht.count()

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals

    # Get the telomere and centromere positions
    telomere_positions, centromere_positions = get_telomere_and_centromere_start_end_positions(tc_intervals_filepath)

    intervals = {}
    for chrom in CHROMS:
        intervals[chrom] = get_intervals_for_chr(
            ht=ht.filter(ht.locus.contig == chrom),
            interval_size=interval_size,
            telomere_positions=telomere_positions[chrom],
            centromere_position=centromere_positions[chrom],
        )

    # Evaluate the intervals
    interval_positions_by_chrom = write_intervals_json(intervals, output_intervals_path)

    return interval_positions_by_chrom


def get_intervals_for_chr(
    ht: hl.Table,
    interval_size: int,
    telomere_positions: dict[str, tuple[int, int]],
    centromere_position: tuple[int, int],
) -> list[hl.Table]:
    """
    Evenly split the chromosome subset of the table into intervals of interval_size.
    Hard exclude the telomeric and centromeric regions of the chromosome.
    """
    # The first valid region for genotyping begins after the first telomere and ends before the centromere
    genotyping_region_1_start = telomere_positions['first'][1] + 1
    genotyping_region_1_end = centromere_position[0] - 1
    ht_region_1 = ht.filter(
        (ht.locus.position >= genotyping_region_1_start) & (ht.locus.position < genotyping_region_1_end),
    )
    ht_region_1 = ht_region_1.add_index(name='row_idx')
    ht_region_1 = ht_region_1.select(ht_region_1.global_row_idx, ht_region_1.row_idx)

    # The second valid region for genotyping starts after the centromere and ends before the second telomere
    genotyping_region_2_start = centromere_position[1] + 1
    genotyping_region_2_end = telomere_positions['second'][0] - 1
    ht_region_2 = ht.filter(
        (ht.locus.position >= genotyping_region_2_start) & (ht.locus.position < genotyping_region_2_end),
    )
    ht_region_2 = ht_region_2.add_index(name='row_idx')
    ht_region_2 = ht_region_2.select(ht_region_2.global_row_idx, ht_region_2.row_idx)

    # Iterate through both genotyping regions in intervals of interval_size, using the row_idx
    intervals = []
    for ht_region in [ht_region_1, ht_region_2]:
        for i in range(0, ht_region.count(), interval_size):
            # Get the interval
            interval = ht_region.filter((ht_region.row_idx >= i) & (ht_region.row_idx < i + interval_size))

            # Add the interval to the intervals list
            intervals.append(interval)

    return intervals


def extract_interval_positions(ht: hl.Table, n_intervals: int) -> List[IntervalInfo]:
    start_time = time.time()
    logging.info(f"Starting interval extraction for {n_intervals} intervals")
   # Get the number of rows in the table
    n_rows = ht.count()

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals

    # Create a new table with interval information
    interval_table = hl.Table.parallelize(
        hl.range(n_intervals).map(lambda i: hl.struct(
            interval_idx = i,
            start_idx = i * interval_size,
            end_idx = hl.min((i + 1) * interval_size - 1, ht.count() - 1),
        )),
    )

    logging.info("Joining interval table with original table")
    # Join the interval table with our original table
    joined = interval_table.key_by(
        global_row_idx = hl.int64(hl.if_else(
            interval_table.interval_idx == n_intervals - 1,
            interval_table.end_idx,
            interval_table.start_idx,
        )),
    ).join(ht.key_by('global_row_idx'))
    logging.info(f"Total rows: {n_rows}, Interval size: {interval_size}")

    # Aggregate to get start and end information for each interval
    # Add a progress column
    joined = joined.annotate(progress = hl.scan.count())

    logging.info("Starting aggregation and collection")
    result = joined.group_by(joined.interval_idx).aggregate(
        start_idx = hl.agg.min(joined.start_idx),
        end_idx = hl.agg.max(joined.end_idx),
        start_pos = hl.agg.min(joined.locus.position),
        end_pos = hl.agg.max(joined.locus.position),
        start_contig = hl.or_missing(hl.agg.count() > 0, hl.agg.take(joined.locus.contig, 1)[0]),
        end_contig = hl.or_missing(hl.agg.count() > 0, hl.agg.take(joined.locus.contig, -1)[0]),
        progress = hl.agg.max(joined.progress),
    )
    # Add action to log progress
    def log_progress(t):
        elapsed_time = time.time() - start_time
        progress_fraction = t.progress / n_rows
        estimated_total_time = elapsed_time / progress_fraction if progress_fraction > 0 else 0
        remaining_time = estimated_total_time - elapsed_time
        logging.info(f"Progress: {progress_fraction:.2%}, Elapsed time: {elapsed_time:.2f}s, Estimated remaining time: {remaining_time:.2f}s")
        return t

    # Collect results with progress logging
    result = result.map(log_progress).collect()

    logging.info("Collection complete. Converting to IntervalInfo objects")
    # Convert to list of IntervalInfo objects
    interval_info = [IntervalInfo(**{k: v for k, v in interval.items() if k != 'progress'}) for interval in result]

    end_time = time.time()
    logging.info(f"Interval extraction completed in {end_time - start_time:.2f} seconds")

    return interval_info


def write_intervals_json(intervals: dict[str, list[hl.Table]], output_path: str):
    """
    Evaluate the list of intervals and save them to a json.
    """
    interval_positions: dict[str, list[tuple[int, int]]] = {}
    for chrom, chrom_intervals in intervals.items():
        interval_positions[chrom] = []
        for i, interval in enumerate(chrom_intervals):
            print(f'Chrom {chrom}, interval {i}: {interval.count_rows()} rows')
            # Get the start and end position of the interval
            positions = interval.aggregate(interval_positions=hl.agg.collect(interval.locus.position))
            start_pos = positions[0]
            end_pos = positions[-1]
            print(f'Interval start: {start_pos}, end: {end_pos}')
            interval_positions[chrom].append((start_pos, end_pos))

    # Save the interval positions to a file
    with open(output_path, 'w') as f:
        json.dump(interval_positions, f)

    print(f'Intervals saved to {output_path}')
    return interval_positions


# def write_intervals_file(output_path: str):
#     """Evaluate the start and end of each interval and write to a file."""
    # Create a temporary Hail Table with a single row
    # temp_ht = hl.Table.parallelize([{'dummy': 0}])
    # temp_ht = temp_ht.annotate(intervals=intervals)
    # # evaluate the intervals
    # result = temp_ht.aggregate(
    #     hl.agg.array_agg(
    #         lambda interval: hl.struct(
    #             start=interval.start,
    #             end=interval.end,
    #         ),
    #         temp_ht.intervals,
    #     ),
    # )
    # with open(output_path, 'w') as f:
    #     for interval in result:
    #         f.write(f'{interval["start"]}\t{interval["end"]}\n')
    #         print(f'{interval["start"]}\t{interval["end"]}')

    # print(f'{len(result)} intervals written to {output_path}')


def main(args):
    """
    Run the script.
    """
    init_batch()
    mt = hl.read_matrix_table(args.input_mt)
    print('Read matrixtable')
    ht = mt.rows()
    ht = ht.select().select_globals()
    ht = ht.add_index(name='global_row_idx')
    ht = ht.key_by('locus', 'alleles', 'global_row_idx')
    print('Re-keyed table')
    with to_path(args.output_intervals_path).open('w') as f:
        for interval in extract_interval_positions(ht, n_intervals=args.n_intervals):
            print(f"Interval: {interval.start_idx}-{interval.end_idx}, "
            f"Positions: {interval.start_pos}-{interval.end_pos}, "
            f"Contigs: {interval.start_contig}-{interval.end_contig}")
            f.write(f"{interval.start_contig}\t{interval.start_pos}\t{interval.end_pos}\t{interval.end_contig}\n")
    # get_even_intervals_from_mt_rows(ht, args.n_intervals, args.tc_intervals_filepath, args.output_intervals_path)
    # intervals = mt._calculate_new_partitions(args.n_intervals)
    # print('Calculated intervals')
    # write_intervals_file(intervals, args.output_intervals_path)


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
