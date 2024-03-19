"""
takes a Pedigree file, generates a Ploidy file from it
"""

from collections import defaultdict

from cpg_utils import to_path


def convert_ped_record(ped_record: str, contigs: list[str], chr_x: str = 'chrX', chr_y: str = 'chrY') -> str:
    """
    Converts a ped file record to a table record.
    Args:
        ped_record (str): line from the PED file
        contigs (str): ordered list of contigs
        chr_x (str): form to use for Chromosome X
        chr_y (str): form to use for Chromosome Y

    Returns:
        ploidy table record
    """

    tokens = ped_record.strip().split('\t')
    sample = tokens[1]
    sex = tokens[4]
    ploidy = defaultdict(lambda: 2)
    if sex == '1':
        ploidy[chr_x] = 1
        ploidy[chr_y] = 1
    elif sex == '2':
        ploidy[chr_x] = 2
        ploidy[chr_y] = 0
    else:
        ploidy[chr_x] = 0
        ploidy[chr_y] = 0
    return '\t'.join([sample] + [str(ploidy[c]) for c in contigs])


def generate_ploidy_table(ped: str, contigs: str, outpath: str):
    """
    Creates a ploidy file from a PED file
    Args:
        ped (str): path to the PED file for this cohort
        contigs (str): path to a file containing the contigs
        outpath (str): path to write the output pploidy table

    Returns:
        None, writes new ploidy table directly
    """

    # read the contigs from the file
    with to_path(contigs).open('r') as f:
        contig_strings = [line.strip().split('\t')[0] for line in f]

    # start the collection of output lines
    output_lines = ['\t'.join(['SAMPLE'] + contig_strings) + '\n']

    # read the PED file
    with to_path(ped).open('r') as ped:
        for line in ped:
            if line.startswith('#'):
                # skip comments / headers
                continue
            output_lines.append(convert_ped_record(ped_record=line, contigs=contig_strings) + '\n')

    # write the data to a file
    with to_path(outpath).open('w') as handle:
        handle.writelines(output_lines)
