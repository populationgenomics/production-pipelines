"""
Read a matrixtable and split it into even intervals of rows.
Use the contigs at the start of each interval to define the 
intervals saved in the output interval_list file.

Hard breaks at centromeres.
"""

import hail as hl
import argparse


def get_even_intervals_from_mt(mt, n_intervals, interval_list_file):
    """
    Get even intervals from a matrixtable.
    """
    # Get the number of rows in the matrixtable
    n_rows = mt.count_rows()

    # Get the number of rows in each interval
    interval_size = n_rows // n_intervals

    # Get the contig and position of the start of each interval
    intervals = []
    for i in range(n_intervals):
        if i == 0:
            start_locus = mt.locus.collect()[0]
        else:
            start_locus = mt.rows().filter(mt.locus.position == 1).collect()[i * interval_size]
        intervals.append(start_locus)

    # Save the intervals to a file
    with open(interval_list_file, 'w') as f:
        for interval in intervals:
            f.write(f"{interval.contig}\t{interval.position}\n")

    return intervals