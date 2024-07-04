"""
Reads all the DeterminePloidy results
Uses the contained ploidy to update the default pedigree
Writes a new pedigree file to use in gCNV
Writes a file containing the SGIDs of samples with aneuploidies
"""

from argparse import ArgumentParser
from os import listdir
from os.path import join


def find_sex(x_ploidy: int | None, y_ploidy: int | None) -> str:
    """
    take the X & Y ploidies, and determine the PED sex
    Args:
        x_ploidy ():
        y_ploidy ():

    Returns:
        int: the value to insert
        1 = male
        2 = female
        0 = unknown/other
    """

    if x_ploidy == 1 and y_ploidy == 1:
        return '1'
    if x_ploidy == 2 and y_ploidy == 0:
        return '2'
    # unacceptableeeee!
    return '0'


def wrangle_genotypes(ploidy_folder: str) -> tuple[dict[str, str], set[str]]:
    """
    trawl the folder for ploidy files, and return a dict

    Args:
        ploidy_folder ():

    Returns:
        a dictionary of SGID: new sex
        a set of all SG IDs where the sample is not cleanly XX or XY
    """

    new_sex_dict: dict[str, str] = {}
    aneuploid_samples: set[str] = set()
    for sample_folder in listdir(ploidy_folder):
        target_file = join(ploidy_folder, sample_folder, 'contig_ploidy.tsv')
        with open(target_file) as handle:
            file_x: int | None = None
            file_y: int | None = None
            topline = handle.readline()
            sgid = topline.rstrip().split('SM:')[1]

            # iterate over everything else
            for line in handle:
                l_list = line.split()
                if l_list[0] == 'chrX':
                    file_x = int(l_list[1])
                if l_list[0] == 'chrY':
                    file_y = int(l_list[1])

            new_sex = find_sex(file_x, file_y)

            if new_sex == '0':
                aneuploid_samples.add(sgid)
            new_sex_dict[sgid] = find_sex(file_x, file_y)

    return new_sex_dict, aneuploid_samples


def update_pedigree(original_ped: str, new_ped: str, new_sexes: dict[str, str]):
    """
    takes the pedigree file, updates sexes, writes back out
    iterate over each line, update the sex col, re-join by tabs

    Args:
        original_ped (str):
        new_ped (str):
        new_sexes (dict[str, int]):
    """

    with open(original_ped) as handle:
        with open(new_ped, 'w') as out_handle:
            for line in handle:
                l_list = line.split()
                sgid = l_list[1]
                l_list[5] = new_sexes[sgid]
                out_handle.write('\t'.join(l_list) + '\n')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('original_ped')
    parser.add_argument('new_ped')
    parser.add_argument('aneuploidies')
    parser.add_argument('ploidy_folder')
    args = parser.parse_args()

    new_sexes, aneuploidies = wrangle_genotypes(args.ploidy_folder)
    print(new_sexes)
    print(aneuploidies)

    update_pedigree(args.original_ped, args.new_ped, new_sexes)

    # open and write, even if the set is empty this will create the file
    with open(args.aneuploidies, 'w') as handle:
        for sample in aneuploidies:
            handle.write(sample + '\n')
