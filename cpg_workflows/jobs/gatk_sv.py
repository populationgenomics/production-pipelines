"""
method(s) for GATK SV PythonJobs
"""


def rename_sv_ids(input_tmp: str, output_file: str):
    """
    A Python method to call as a PythonJob, edits content of the VCF
    Replaces the standard ID with a guaranteed-unique ID (hopefully)
    The new ID compounds Type, Chr & Start, then either end or Chr2 & End2

    Writes the file back out to the specified path

    Args:
        input_tmp (str): path to temp file generated by merging
        output_file (str): path to write uncompressed edited version to
    """
    import gzip

    headers = []
    others = []

    # crack open that VCF and have a little sip
    with gzip.open(input_tmp, 'rt') as f:
        for line in f:
            # don't alter current header lines
            if line.startswith('#'):
                headers.append(line)
                continue

            # split on tabs
            l_split = line.split('\t')

            # e.g. AN_Orig=61;END=56855888;SVTYPE=DUP
            info_dict: dict[str, str] = dict(el.split('=') for el in l_split[7].split(';'))

            # strip out the "chr" prefix to abbreviate String
            chrom = l_split[0].removeprefix('chr')
            start = l_split[1]

            # e.g. <DEL> -> DEL
            alt_allele = l_split[4][1:-1]

            # maybe there's a second chromosome?
            if all(key in info_dict for key in ['CHR2', 'END2']):
                chrom2 = info_dict['CHR2'].removeprefix('chr')
                end2 = info_dict['END2']
                # e.g. BND_1-12345_2-67890
                l_split[2] = f'{alt_allele}_{chrom}-{start}_{chrom2}-{end2}'

            else:
                # e.g. CNV_1-12345-67890
                l_split[2] = f'{alt_allele}_{chrom}-{start}-{info_dict["END"]}'

            # rebuild the line
            others.append('\t'.join(l_split))

    with open(output_file, 'w') as f:
        f.writelines(headers)
        f.writelines(others)