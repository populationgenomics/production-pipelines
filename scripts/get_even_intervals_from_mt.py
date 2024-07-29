"""
Read a matrixtable and split it into even intervals of rows.
Use the contigs at the start of each interval to define the 
intervals saved in the output interval_list file.

Hard breaks at centromeres.
"""

import hail as hl
import argparse


def get_intervals_for_chr()

def get_even_intervals_from_mt(mt: hl.MatrixTable, n_intervals: int, interval_list_file):
    """
    Get even intervals from a matrixtable.
    """
    # Get the number of rows in the matrixtable        
    n_rows = mt.count_rows()

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals

    for chrom in ['chr1',]:
        get_intervals_for_chr(mt, chrom, n_intervals, interval_size, centromere_position, interval_list_file)

    # Add the integer index of each row as a new row field.
    mt = mt.add_row_index(name='row_idx')
    mt = mt.rows().select('row_idx', 'locus')

    # Iterate through the mt in intervals of interval_size, using the row_idx
    intervals = []
    for n in range(n_intervals):
        # Get the start and end of the interval
        start = n * interval_size
        # stop = (n + 1) * interval_size - 1

        intervals.append(
            mt.filter_rows(mt.row_idx == start)['locus'],
        )

    # Evaluate the contig and position of each interval 
    intervals = [interval.collect()[0] for interval in intervals]

    # Save the intervals to a file
    # with open(interval_list_file, 'w') as f:
    #     for interval in intervals:
    #         f.write(f"{interval.contig}\t{interval.position}\n")

    return intervals