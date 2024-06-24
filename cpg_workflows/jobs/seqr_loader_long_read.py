def modify_sniffles_vcf(file_in: str, file_out: str):
    """
    to be run as a PythonJob - scrolls through the

    Args:
        file_in (str): localised, VCF directly from Sniffles
        file_out (str): local batch output path, same VCF with INFO/ALT alterations
    """
    import gzip

    # read and write compressed. This is only a single sample VCF, but... it's good practice
    with gzip.open(file_in, 'rt') as f:
        with gzip.open(file_out, 'wt') as f_out:
            for line in f:
                # don't alter current header lines
                if line.startswith('#'):
                    f_out.write(line)
                    continue

                # for non-header lines, split on tabs
                l_split = line.split('\t')

                # reduce the massive REF alleles to a single base
                l_split[3] = l_split[3][0]

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
