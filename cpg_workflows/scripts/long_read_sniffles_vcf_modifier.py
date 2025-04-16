import gzip
from argparse import ArgumentParser

from pyfaidx import Fasta

DUP = 'DUP'
DEL = 'DEL'
CHRX = 'chrX'
CHRY = 'chrY'
CHRM = 'chrM'

SV_ANNOTATION_TYPES = [
    # From GATKSVVCFConstants.StructuralVariantAnnotationType enum
    'BND',
    'CNV',
    'CPX',
    'CTX',
    'DEL',
    'DUP',
    'INS',
    'INV',
]


def read_sex_mapping_file(sex_mapping_file_path: str) -> dict[str, int]:
    """
    Read in the mapping file and return the mapping of LRS ID to sex value
    """
    sex_mapping = {}
    with open(sex_mapping_file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            sex_mapping[line.split()[0]] = int(line.split()[1])

    return sex_mapping


def translate_var_and_sex_to_cn(contig: str, var_type: str, genotype: str, sex: int) -> int | str:
    """
    Translate a variant and sex to a CN value
    using CN==2 as a baseline, this is modified up or down based on the variant call
    pull out the numbers from the genotype String
    only need a CN for Deletion and Duplication

    Args:
        contig ():
        var_type ():
        genotype (): GT string, e.g. 0/1, 0|1, 1|0
        sex (int): 0=Unknown, 1=Male, 2=Female

    Returns:
        int | string, the CN value (copy number) or '.' if not applicable
    """
    if '0' not in genotype and '1' not in genotype:
        # if there are no 0s or 1s, we can't determine the CN
        return '.'

    # determine the baseline copy number
    # conditions where the copy number is 1
    if sex == 1 and contig in [CHRX, CHRY]:
        copy_number = 1
    elif sex == 2 and contig == CHRY:
        copy_number = 0
    else:
        copy_number = 2

    if (var_type not in (DUP, DEL)) or contig == CHRM:
        return copy_number

    # find the number of copies of the variant. Normalised variants, so only one alt
    num_alt = genotype.count('1')

    # reduce or increase the copy number based on the genotype, depending on the variant type
    if var_type == DUP:
        copy_number += num_alt
    elif var_type == DEL:
        copy_number -= num_alt
    return copy_number


def cli_main():
    parser = ArgumentParser(description='CLI for the Sniffles VCF modification script')
    parser.add_argument('--vcf_in', help='Path to a localised VCF, this will be modified', required=True)
    parser.add_argument('--vcf_out', help='Path to an output location for the modified VCF', required=True)
    parser.add_argument('--fa', help='Path to a FASTA sequence file for GRCh38', required=True)
    parser.add_argument(
        '--sex_mapping_file',
        help='Path to a file containing LRS IDs and participant sexes',
        required=True,
    )
    args = parser.parse_args()
    modify_sniffles_vcf(
        file_in=args.vcf_in,
        file_out=args.vcf_out,
        fa=args.fa,
        sex_mapping_file=args.sex_mapping_file,
    )


def modify_sniffles_vcf(file_in: str, file_out: str, fa: str, sex_mapping_file: str):
    """
    Scrolls through the SV VCF and performs a few updates:
    - replaces the ALT allele with a symbolic "<TYPE>", derived from the SVTYPE INFO field
    - swaps out the REF (huge for deletions, a symbolic "N" for insertions) with the ref base
    - adds a CN field to the FORMAT field, and fills it with copy number values based on the genotype/sex

    rebuilds the VCF following those edits, and writes the compressed data back out

    Args:
        file_in (str): localised, VCF directly from Sniffles
        file_out (str): local batch output path, same VCF with INFO/ALT alterations
        fa (str): path to a reference FastA file, requires an implicit fa.fai index
        sex_mapping_file (str): path to a file which maps LRS ID to sex value
    """

    # as_raw as a specifier here means that get_seq queries are just the sequence, no contig ID attached
    fasta_client = Fasta(filename=fa, as_raw=True)

    # Get the dict of LRS ID to sex value
    sex_mapping = read_sex_mapping_file(sex_mapping_file)

    # read and write compressed. This is only a single sample VCF, but... it's good practice
    with gzip.open(file_in, 'rt') as f, gzip.open(file_out, 'wt') as f_out:

        for index, line in enumerate(f):

            if index % 10000 == 0:
                print(f'Lines processed: {index}')

            # alter the sample line in the header
            if line.startswith('#'):
                if 'FORMAT=<ID=GT' in line:
                    f_out.write('##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description="Copy number of this variant">\n')

                if 'CHROM' in line:
                    # Find the sample IDs in the header line to match to the input sex values
                    l_split = line.rstrip().split('\t')
                    sample_ids = l_split[9:]

                f_out.write(line)
                continue

            # for non-header lines, split on tabs
            l_split = line.rstrip().split('\t')

            # set the reference allele to be the correct reference base
            chrom = l_split[0]
            position = int(l_split[1])
            try:
                new_base = fasta_client.get_seq(chrom, position, position)
            except Exception as e:
                print(f'Error getting sequence for {chrom}:{position} - {e}')
                continue

            # a quick check, if we can
            if l_split[3] not in ('N', '.'):
                # If using the hg38 masked reference, the base will always be upper case
                # So make sure the comparison to the Sniffles REF is upper case too
                assert (
                    new_base == l_split[3][0].upper()
                ), f'Discrepancy between faidx and Sniffles: {new_base}, {l_split[3]}'

            # replace the REF String
            l_split[3] = new_base

            # e.g. AN_Orig=61;END=56855888;SVTYPE=DUP
            info_dict: dict[str, str] = {}
            for entry in l_split[7].split(';'):
                if '=' in entry:
                    key, value = entry.split('=')
                    info_dict[key] = value

            # get the SVTYPE, always present
            sv_type = info_dict['SVTYPE']
            if sv_type not in SV_ANNOTATION_TYPES:
                print(f'Forbidden SVTYPE: {sv_type}, at {chrom}:{position}')
                continue

            # update the FORMAT schema field
            l_split[8] = f'{l_split[8]}:RD_CN'

            for i, sample_id in enumerate(sample_ids, start=9):
                # pull out the GT section of the FORMAT field
                gt_string = dict(zip(l_split[8].split(':'), l_split[i].split(':')))['GT']

                # determine the copy number, based on deviation from the baseline
                copy_number = translate_var_and_sex_to_cn(
                    contig=chrom,
                    var_type=sv_type,
                    genotype=gt_string,
                    sex=sex_mapping[sample_id],
                )

                # update the FORMAT content field for this sample
                l_split[i] = f'{l_split[i]}:{copy_number}'

            # replace the alt with a symbolic String
            l_split[4] = f'<{sv_type}>'

            # breakends (BND) aren't annotated with an END or SVLEN, so we use the CHR2 value
            if 'END' in info_dict:
                end_position = info_dict['END']
            elif 'CHR2' in info_dict:
                # No END in INFO, using CHR2
                end_position = info_dict['CHR2']
            else:
                end_position = str(position)

            # replace the UID with something meaningful: type_chrom_pos_end
            # this is required as GATK-SV's annotation module sorts on ID, not on anything useful
            l_split[2] = f'{sv_type}_{chrom}_{position}_{end_position}'

            # rebuild the string and write as output
            f_out.write('\t'.join(l_split) + '\n')


if __name__ == '__main__':
    # if called as a script, call through to the ArgParse CLI
    cli_main()
