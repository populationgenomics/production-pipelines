import gzip
from argparse import ArgumentParser

from pyfaidx import Fasta


def cli_main():
    parser = ArgumentParser(description='CLI for the Sniffles VCF modification script')
    parser.add_argument('--vcf_in', help='Path to a localised VCF, this will be modified', required=True)
    parser.add_argument('--vcf_out', help='Path to an output location for the modified VCF', required=True)
    parser.add_argument('--fa', help='Path to a FASTA sequence file for GRCh38', required=True)
    parser.add_argument('--fai', help='Path to a FASTA.fai sequence index file for GRCh38', required=True)
    parser.add_argument('--ext_id', help='Path to the Sample ID in the input VCF', default=None)
    parser.add_argument('--int_id', help='Path to the Sample ID we want in the output VCF', default=None)
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(f'Unknown input flag(s) used: {unknown}')

    modify_sniffles_vcf(
        file_in=args.vcf_in,
        file_out=args.vcf_out,
        fa=args.fa,
        fa_fai=args.fai,
        ext_id=args.ext_id,
        int_id=args.int_id,
    )


def modify_sniffles_vcf(
    file_in: str,
    file_out: str,
    fa: str,
    fa_fai: str,
    ext_id: str | None = None,
    int_id: str | None = None,
):
    """
    Scrolls through the VCF and performs a few updates:

    - replaces the External Sample ID with the internal CPG identifier (if both are provided)
    - replaces the ALT allele with a symbolic "<TYPE>", derived from the SVTYPE INFO field
    - swaps out the REF (huge for deletions, a symbolic "N" for insertions) with the ref base

    rebuilds the VCF following those edits, and writes the compressed data back out

    Args:
        file_in (str): localised, VCF directly from Sniffles
        file_out (str): local batch output path, same VCF with INFO/ALT alterations
        fa (str): path to a reference FastA file
        fa_fai (str): path to the FA index
        ext_id (str): external ID to replace (if found)
        int_id (str): CPG ID, required inside the reformatted VCF
    """

    fasta_client = Fasta(filename=fa, indexname=fa_fai, as_raw=True)

    # read and write compressed. This is only a single sample VCF, but... it's good practice
    with gzip.open(file_in, 'rt') as f, gzip.open(file_out, 'wt') as f_out:

        for index, line in enumerate(f):

            if index % 100 == 0:
                print(f'line {index}')
            # alter the sample line in the header
            if line.startswith('#'):
                if line.startswith('#CHR') and (ext_id and int_id):
                    line.replace(ext_id, int_id)

                f_out.write(line)
                continue

            # for non-header lines, split on tabs
            l_split = line.split('\t')

            # set the reference allele to be the correct reference base
            new_base = fasta_client.get_seq(l_split[0], int(l_split[1]), int(l_split[1]))

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

            # replace the alt with a symbolic String
            l_split[4] = f'<{info_dict["SVTYPE"]}>'
            # rebuild the string and write as output
            f_out.write('\t'.join(l_split))


if __name__ == '__main__':
    # if called as a script, call through to the ArgParse CLI
    cli_main()
