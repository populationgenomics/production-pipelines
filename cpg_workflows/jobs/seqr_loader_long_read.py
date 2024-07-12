def modify_sniffles_vcf(file_in: str, file_out: str, ext_id: str, int_id: str, fa: str, fa_fai: str):
    """
    to be run as a PythonJob - scrolls through the VCF and performs a few updates:

    - replaces the External Sample ID with the internal CPG identifier
    - replaces the ALT allele with a symbolic "<TYPE>", derived from the SVTYPE INFO field
    - swaps out the REF (huge for deletions, a symbolic "N" for insertions) with the ref base

    rebuilds the VCF following those edits, and writes the compressed data back out

    Args:
        file_in (str): localised, VCF directly from Sniffles
        file_out (str): local batch output path, same VCF with INFO/ALT alterations
        ext_id (str): external ID to replace (if found)
        int_id (str): CPG ID, required inside the reformatted VCF
        fa (str): path to a reference FastA file
        fa_fai (str): path to the index
    """
    import gzip

    import hail as hl

    # initiate a batch
    hl.init()
    # set the default reference
    hl.default_reference('GRCh38')

    # create a ReferenceGenome object
    rg_38 = hl.ReferenceGenome('GRCh38')

    # add the sequence
    rg_38.add_sequence(fasta_file=fa, index_file=fa_fai)

    # read and write compressed. This is only a single sample VCF, but... it's good practice
    with gzip.open(file_in, 'rt') as f, gzip.open(file_out, 'wt') as f_out:

        for line in f:
            # alter the sample line in the header
            if line.startswith('#'):
                if line.startswith('#CHR'):
                    line.replace(ext_id, int_id)

                f_out.write(line)
                continue

            # for non-header lines, split on tabs
            l_split = line.split('\t')

            # set the reference allele to be the correct reference base
            l_split[3] = hl.eval(hl.get_sequence(l_split[0], int(l_split[1]), reference_genome=rg_38))

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
