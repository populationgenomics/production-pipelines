import logging

import hail as hl

from cpg_utils.config import reference_path

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_cadd_grch38_ht(outpath: str, tmp_dir: str) -> hl.Table:
    """
    Create a Hail Table with CADD scores for GRCh38.

    The combined CADD scores in the returned table are from the following sources:
        - all SNVs: `cadd.v1.6.whole_genome_SNVs.tsv.bgz` (81G) downloaded from
          `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz`.
          It contains 8,812,917,339 SNVs.
        - gnomad 3.0 indels: `cadd.v1.6.gnomad.genomes.v3.0.indel.tsv.bgz` (1.1G)
          downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`.
          It contains 100,546,109 indels from gnomaD v3.0.
        - gnomad 3.1 indels: `cadd.v1.6.gnomad.genomes.v3.1.indels.new.ht` was run
          on gnomAD v3.1 with CADD v1.6 in 2020. It contains 166,122,720 new indels from
          gnomAD v3.1 compared to v3.0.
        - gnomad 3.1 complex indels: `cadd.v1.6.gnomad.genomes.v3.1.indels.complex.ht`
          was run on gnomAD v3.1 with CADD v1.6 in 2020. It contains 2,307 complex
          variants that do not fit Hail's criteria for an indel and thus exist in
          a separate table than the gnomad 3.1 indels.
        - gnomAD v4 exomes indels: `cadd.v1.6.gnomad.exomes.v4.0.indels.new.tsv.bgz`
          (368M) was run on gnomAD v4 with CADD v1.6 in 2023. It contains 32,561,
          253 indels that are new in gnomAD v4.
        - gnomAD v4 genomes indels: `cadd.v1.6.gnomad.genomes.v4.0.indels.new.tsv.bgz`
          (13M) was run on gnomAD v4 with CADD v1.6 in 2023. It contains 904,906 indels
          that are new in gnomAD v4 genomes because of the addition of HGDP/TGP samples.

    .. note::

         ~1,9M indels were duplicated in gnomAD v3.0 and v4.0 or in gnomAD v3.1 and
         v4.0. However, CADD only generates a score per loci. We keep only the latest
         prediction, v4.0, for these loci.
         The output generated a CADD HT with 9,110,177,520 rows.

    :return: Hail Table with CADD scores for GRCh38.
    """

    def _load_cadd_raw(cadd_tsv, tmp_dir: str, tmp_file: str) -> hl.Table:
        """Functions to load CADD raw data in TSV format to Hail Table."""
        column_names = {
            "f0": "chr",
            "f1": "pos",
            "f2": "ref",
            "f3": "alt",
            "f4": "RawScore",
            "f5": "PHRED",
        }
        types = {"f0": hl.tstr, "f1": hl.tint32, "f4": hl.tfloat32, "f5": hl.tfloat32}
        ht = hl.import_table(
            cadd_tsv,
            types=types,
            no_header=True,
            force_bgz=True,
            comment="#",
            min_partitions=1000,
        )
        ht = ht.rename(column_names)
        chr = hl.if_else(ht.chr.startswith("chr"), ht.chr, hl.format("chr%s", ht.chr))
        ht = ht.annotate(
            locus=hl.locus(chr, ht.pos, reference_genome="GRCh38"),
            alleles=hl.array([ht.ref, ht.alt]),
        )
        ht = ht.select("locus", "alleles", "RawScore", "PHRED")
        ht = ht.key_by("locus", "alleles")
        # ht = ht.checkpoint(new_temp_file("cadd", "ht"))
        logger.info(f'Checkpointing to {tmp_dir}/{tmp_file}')
        return ht.checkpoint(tmp_dir + '/' + tmp_file, overwrite=True)

    logging.info(f'CADD path: {reference_path("CADD_v1.7_snvs")}')
    snvs = _load_cadd_raw(
        reference_path('CADD_v1.7_snvs'),
        tmp_dir,
        'snvs.ht',
    )

    logging.info(f'gnomad v3.0 indels path: {reference_path("exomiser_cadd/indel_tsv")}')
    indel3_0 = _load_cadd_raw(
        reference_path('exomiser_cadd/indel_tsv'),
        tmp_dir,
        'indels3_0.ht',
    )

    # TODO: Not sure if these are actually released, so not using them for now.
    # indel3_1 = hl.read_table(
    #     reference_path('gnomad.genomes.v3.1.indels.new.ht'),
    # )
    # indel3_1_complex = hl.read_table(
    #     reference_path('gnomad.genomes.v3.1.indels.complex.ht'),
    # )

    # TODO: upload v4.0 exome indels
    # indel4_e = _load_cadd_raw(
    #     reference_path('gnomad.exomes.v4.0.indels.new.tsv.bgz'),
    # )
    logging.info(f'gnomad v4.0 genomes indels path: {reference_path("CADD_v1.7_indels")}')
    indel4_g = _load_cadd_raw(
        reference_path('CADD_v1.7_indels'),
        tmp_dir,
        'indels4_g.ht',
    )

    # Merge the CADD predictions run for v3 versions.
    # indel3 = indel3_0.union(indel3_1, indel3_1_complex)

    # Merge the CADD predictions run for v4 versions.
    # indel4 = indel4_e.union(indel4_g).distinct()

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3_0.anti_join(indel4_g)

    ht = snvs.union(indel3, indel4_g)

    ht = ht.select(cadd=hl.struct(phred=ht.PHRED, raw_score=ht.RawScore))
    ht = ht.annotate_globals(cadd_version="v1.7")

    ht = ht.checkpoint(outpath, overwrite=True)
    return ht
