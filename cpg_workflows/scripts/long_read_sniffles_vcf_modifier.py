import gzip
from argparse import ArgumentParser

from pyfaidx import Fasta

DUP = 'DUP'
DEL = 'DEL'
CHRX = 'chrX'
CHRY = 'chrY'
CHRM = 'chrM'


def translate_var_and_sex_to_cn(contig: str, var_type: str, genotype: str, sex: int) -> int:
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
        int, the CN value (copy number)
    """

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
    parser.add_argument('--new_id', help='The Sample ID we want in the output VCF', default=None)
    parser.add_argument('--fa', help='Path to a FASTA sequence file for GRCh38', required=True)
    parser.add_argument('--sex', help='0=Unknown,1=Male, 2=Female', default=0, type=int)
    parser.add_argument(
        '--sv',
        help='Boolean flag to indicate if the VCF is SVs. False=SNPs_Indels',
        action='store_true',
    )
    args = parser.parse_args()

    modify_sniffles_vcf(
        file_in=args.vcf_in,
        file_out=args.vcf_out,
        fa=args.fa,
        new_id=args.new_id,
        sex=args.sex,
        sv=args.sv,
    )


def modify_sniffles_vcf(file_in: str, file_out: str, fa: str, new_id: str | None = None, sex: int = 0, sv: bool = True):
    """
    Scrolls through the VCF and performs a few updates:

    - replaces the External Sample ID with the internal CPG identifier (if both are provided)

    And if the VCF contains SVs:
    - replaces the ALT allele with a symbolic "<TYPE>", derived from the SVTYPE INFO field
    - swaps out the REF (huge for deletions, a symbolic "N" for insertions) with the ref base

    rebuilds the VCF following those edits, and writes the compressed data back out

    Args:
        file_in (str): localised, VCF directly from Sniffles
        file_out (str): local batch output path, same VCF with INFO/ALT alterations
        fa (str): path to a reference FastA file, requires an implicit fa.fai index
        new_id (str): CPG ID, required inside the reformatted VCF
        sex (int): 0=Unknown, 1=Male, 2=Female
        sv (bool): True=SV, False=SNPs_Indels
    """

    # as_raw as a specifier here means that get_seq queries are just the sequence, no contig ID attached
    fasta_client = Fasta(filename=fa, as_raw=True)

    # read and write compressed. This is only a single sample VCF, but... it's good practice
    with gzip.open(file_in, 'rt') as f, gzip.open(file_out, 'wt') as f_out:

        for index, line in enumerate(f):

            if index % 10000 == 0:
                print(f'Lines processed: {index}')

            # alter the sample line in the header
            if line.startswith('#'):
                if 'FORMAT=<ID=ID' in line and sv:
                    f_out.write('##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description="Copy number of this variant">\n')

                if line.startswith('#CHR') and new_id:
                    l_split = line.rstrip().split('\t')
                    l_split[9] = new_id
                    f_out.write('\t'.join(l_split) + '\n')
                    continue

                if 'FORMAT=<ID=AF,Number=1' in line and not sv:
                    # Correct the AF field to have 'Number=A' for SNPs/Indels, to allow for multiple AF values
                    line = line.replace('Number=1', 'Number=A')

                f_out.write(line)
                continue

            if not sv:
                # if we're not dealing with SVs, just write the line and move on
                f_out.write(line)
                continue

            # for non-header lines, split on tabs
            l_split = line.rstrip().split('\t')

            # set the reference allele to be the correct reference base
            chrom = l_split[0]
            position = int(l_split[1])
            new_base = fasta_client.get_seq(chrom, position, position)

            # a quick check, if we can
            if l_split[3] != 'N':
                assert new_base == l_split[3][0], f'Discrepancy between faidx and Sniffles: {new_base}, {l_split[3]}'

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

            # pull out the GT section of the FORMAT field
            gt_string = dict(zip(l_split[8].split(':'), l_split[9].split(':')))['GT']

            # determine the copy number, based on deviation from the baseline
            copy_number = translate_var_and_sex_to_cn(contig=chrom, var_type=sv_type, genotype=gt_string, sex=sex)

            # update the FORMAT schema field
            l_split[8] = f'{l_split[8]}:CN'

            # update the FORMAT content field
            l_split[9] = f'{l_split[9]}:{copy_number}'

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
