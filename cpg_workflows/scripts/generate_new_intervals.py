"""
Generate new intervals from a MatrixTable based on data distribution

First we read the MT, and generate rough intervals based on the rows
My understanding is that this is a based on the distribution of keys, so
each partition will have roughly the same number of keys (i.e. variants)

This initial group of intervals is split such that each interval fully sits
on a contig. This is to avoid issues with cross-contig intervals.

Then we polish the intervals by dodging centromeres and telomeres. This is
done by reading in an interval_list file, which contains the centromere and
telomere regions. We then check if the intervals overlap with these regions,
and if so, we split the intervals around the centromere/telomere.

For telomere overlaps, we shift the start or end position of the interval
to outside the telomere region.

For centromere overlaps, we split the intervals around the centromere, creating
two new intervals.

The input argument --intervals is the initial number of intervals created. Due to
splitting this can then jump up a bit (in a small test this was +50). This may be
negligible on larger splits, but for smaller interval counts can be huge.
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch


def get_naive_intervals(mt: hl.MatrixTable, intervals: int) -> list[tuple[str, int, int]]:
    """
    Get naive new intervals from a MatrixTable

    Args:
        mt (hl.MatrixTable): Input MatrixTable
        intervals (int): Number of intervals to generate

    Returns:
        list[tuple[str, int, int]]: List of intervals, in non-hail form
    """
    naive_interval_structs = mt._calculate_new_partitions(intervals)
    logging.info(f'Generated {len(naive_interval_structs)} intervals')

    # store both start and end positions separately
    new_intervals = [
        (
            interval.start.locus.contig,
            interval.start.locus.position,
            interval.end.locus.contig,
            interval.end.locus.position,
        )
        for interval in naive_interval_structs
    ]

    # check for cross-contig intervals, and split
    final_intervals: list[tuple[str, int, int]] = []
    parsed_chroms: set[str] = set()
    for index, (contig1, start, contig2, end) in enumerate(new_intervals, start=1):
        # shift start position to 1, regardless of the variant positions used when generating the intervals
        if contig1 not in parsed_chroms:
            parsed_chroms.add(contig1)
            start = 1

        if contig1 != contig2:
            parsed_chroms.add(contig2)
            # get the end of the first
            first_end = hl.get_reference('GRCh38').lengths[contig1]
            final_intervals.append((contig1, start, first_end))
            final_intervals.append((contig2, 1, end))

        else:
            # if this is the final interval on this chromosome, shift the end position to the end of the chromosome
            # true if this is the latest interval on this chromosome
            if index == len(new_intervals) or new_intervals[index][0] != contig1:
                final_intervals.append((contig1, start, hl.get_reference('GRCh38').lengths[contig1]))
            else:
                final_intervals.append((contig1, start, end))

    logging.info(f'Naive intervals post splitting: {len(final_intervals)}')

    return final_intervals


def overlaps(a: tuple[int, int], b: tuple[int, int]) -> int:
    """
    Return the amount of overlap, in bp
    between a and b.
    If >0, the number of bp of overlap
    If 0,  they are book-ended.
    If <0, the distance in bp between them
    """

    if (overlap := min(a[1], b[1]) - max(a[0], b[0])) > 0:
        return overlap
    return 0


def split_interval(interval: tuple[str, int, int], max_interval_size: int) -> list[tuple[str, int, int]]:
    """
    Split an interval into approximately evenly sized smaller intervals of at most max_interval_size bp
    """
    chrom, start, end = interval
    interval_size = end - start

    num_splits = (interval_size + max_interval_size - 1) // max_interval_size
    split_size = interval_size // num_splits

    new_intervals = []
    for i in range(num_splits):
        new_start = start + i * split_size
        new_end = new_start + split_size if i < num_splits - 1 else end
        new_intervals.append((chrom, new_start, new_end))

    return new_intervals


def polish_intervals(
    naive_intervals: list[tuple[str, int, int]],
    meres_file: str,
    max_length: int = 3000000,
) -> list[tuple[str, int, int]]:
    """
    Polish intervals by splitting/removing centromere and telomere regions
    and setting a max interval length
    Args:
        naive_intervals (list): naive intervals, each limited to a single contig
        meres_file (str): path to file containing centromere and telomere regions

    Returns:
        a polished version of the intervals file, where no centromere or telomere regions are present
    """
    telomeres: dict[str, list[tuple[int, int]]] = {}
    centromeres: dict[str, tuple[int, int]] = {}

    logging.info(f'Reading centromere and telomere regions from {meres_file}')

    with open(meres_file, encoding='utf-8') as f:
        for line in f:
            if line.startswith('@SQ') or line.startswith('@HD'):
                continue

            contig, start_str, end_str, strand, region_type = line.strip().split('\t')
            if 'centromere' in region_type:
                centromeres[contig] = (int(start_str), int(end_str))
            else:
                telomeres.setdefault(contig, []).append((int(start_str), int(end_str)))

    new_intervals: list[tuple[str, int, int]] = []
    logging.info('Polishing intervals')
    for chrom, start, end in naive_intervals:
        found_overlap = False
        # does this overlap with a telomere?
        for telo in telomeres.get(chrom, []):
            # there's an overlap
            if overlaps(telo, (start, end)):
                # this is a 'start' telomere, shift the start coord
                if telo[0] == 1:
                    start = telo[1]
                else:
                    end = telo[0]
        # does it overlap with a centromere? If so split around it
        if centro_region := centromeres.get(chrom):
            # check for an overlap
            if overlaps(centro_region, (start, end)):
                # check we have some interval left
                if centro_region[0] > start:
                    new_intervals.append((chrom, start, centro_region[0]))
                if end > centro_region[1]:
                    new_intervals.append((chrom, centro_region[1], end))
                continue

        # check the length of the interval
        if end - start > max_length:
            # split the interval
            new_intervals.extend(split_interval((chrom, start, end), max_length))
            continue
        new_intervals.append((chrom, start, end))

    logging.info(f'Final intervals: {len(new_intervals)}')

    return new_intervals


def cli_main():
    parser = ArgumentParser(description='Generate new intervals from a list of intervals')
    parser.add_argument('--mt', help='Input MatrixTable')
    parser.add_argument('--out', help='Output intervals file')
    parser.add_argument('--intervals', help='Intervals to generate', type=int, default=100)
    parser.add_argument('--max_length', help='Max length of an interval', type=int, default=3000000)
    parser.add_argument('--meres_file', help='interval_list file of Centromere & telomere regions')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    main(args.mt, args.out, intervals=args.intervals, max_length=args.max_length, meres=args.meres_file)


def main(mt: str, out: str, intervals: int, max_length: int, meres: str):
    init_batch()
    mt = hl.read_matrix_table(mt)

    # generate rough intervals based on the rows in this MT
    new_intervals = get_naive_intervals(mt, intervals)

    # polish the intervals by dodging centromeres and telomeres and capping the length
    better_intervals = polish_intervals(new_intervals, meres, max_length)
    with open(out, 'w') as f:
        for contig, start, end in better_intervals:
            f.write(f'{contig}\t{start}\t{end}\n')


if __name__ == "__main__":
    cli_main()
