import logging

import hail as hl

from cpg_utils.config import reference_path

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

NO_CHR_TO_CHR_CONTIG_RECODING = {
    "1": "chr1",
    "2": "chr2",
    "3": "chr3",
    "4": "chr4",
    "5": "chr5",
    "6": "chr6",
    "7": "chr7",
    "8": "chr8",
    "9": "chr9",
    "10": "chr10",
    "11": "chr11",
    "12": "chr12",
    "13": "chr13",
    "14": "chr14",
    "15": "chr15",
    "16": "chr16",
    "17": "chr17",
    "18": "chr18",
    "19": "chr19",
    "20": "chr20",
    "21": "chr21",
    "22": "chr22",
    "X": "chrX",
    "Y": "chrY",
    "MT": "chrM",
}


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
        logger.info(f'Checkpointing to {tmp_dir}/{tmp_file}')
        return ht.checkpoint(tmp_dir + '/' + tmp_file, overwrite=True)

    snvs = _load_cadd_raw(
        reference_path('CADD_v1.7_snvs'),
        tmp_dir,
        'snvs.ht',
    )

    indel3_0 = _load_cadd_raw(
        reference_path('exomiser_cadd/indel_tsv'),
        tmp_dir,
        'indels3_0.ht',
    )

    # TODO: Check if gnomad.genomes.v3.1.indels.new.ht and gnomad.genomes.v3.1.indels.complex.ht are released and available for use.
    # TODO: Upload and use v4.0
    logger.info(f'gnomad v4.0 genomes indels path: {reference_path("CADD_v1.7_indels")}')
    indel4_g = _load_cadd_raw(
        reference_path('CADD_v1.7_indels'),
        tmp_dir,
        'indels4_g.ht',
    )

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3_0.anti_join(indel4_g)

    ht = snvs.union(indel3, indel4_g)

    ht = ht.select(cadd=hl.struct(phred=ht.PHRED, raw_score=ht.RawScore))
    ht = ht.annotate_globals(cadd_version="v1.7")

    ht = ht.checkpoint(outpath, overwrite=True)
    return ht


def create_spliceai_grch38_ht(output_path: str) -> hl.Table:
    """
    Create a Hail Table with SpliceAI scores for GRCh38.

    SpliceAI scores are from the following resources:
        - Precomputed SNVs: spliceai_scores.masked.snv.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - Precomputed indels: spliceai_scores.masked.indel.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - gnomAD v3 indels: gnomad_v3_indel.spliceai_masked.vcf.bgz,
          computed on v3.1 indels by Illumina in 2020.
        - gnomAD v4 indels: gnomad_v4_new_indels.spliceai_masked.vcf.bgz,
          computed on v4 indels that are new compared to v3 indels by Illumina in
          February 2023.
        - gnomAD v3 and v4 unscored indels:
          spliceai_scores.masked.gnomad_v3_v4_unscored_indels.hg38.vcf.bgz,
          another set of indels were not scored in v3 or v4 but computed by Illumina in
          September 2023.

    :return: Hail Table with SpliceAI scores for GRCh38.
    """
    snvs_path = reference_path('ourdna_browser/splice_ai_snvs')
    indels_path = reference_path('ourdna_browser/splice_ai_indels')

    def import_spliceai_vcf(path: str) -> hl.Table:
        """
        Import SpliceAI VCF into Hail Table.

        :param str path: Path to SpliceAI VCF.
        :return: Hail MatrixTable with SpliceAI scores.
        """
        ht = hl.import_vcf(
            path,
            force_bgz=True,
            # header_file=header_file_path,
            reference_genome="GRCh38",
            contig_recoding=NO_CHR_TO_CHR_CONTIG_RECODING,
            skip_invalid_loci=True,
            min_partitions=1000,
        ).rows()
        return ht

    logger.info("Importing vcf of SpliceAI scores into HT...")
    spliceai_snvs = import_spliceai_vcf(snvs_path)
    spliceai_indels = import_spliceai_vcf(indels_path)

    ht = spliceai_snvs.union(spliceai_indels)

    # `explode` because some variants fall on multiple genes and have a score per gene.
    # All rows without a SpliceAI score, an empty array, are removed through explode.
    logger.info("Exploding SpliceAI scores...")
    ht = ht.explode(ht.info.SpliceAI)

    # delta_score array for 4 splicing consequences: DS_AG|DS_AL|DS_DG|DS_DL
    logger.info("Annotating SpliceAI scores...")
    delta_scores = ht.info.SpliceAI.split(delim="\\|")[2:6]
    logger.info(f"Delta scores: {delta_scores}")
    ht = ht.annotate(delta_scores=hl.map(lambda x: hl.float32(x), delta_scores))

    logger.info(
        "Getting the max SpliceAI score across consequences for each variant per" " gene...",
    )
    ht = ht.select(ds_max=hl.max(ht.delta_scores))

    logger.info("Getting the max SpliceAI score for each variant across genes...")
    ht = ht.collect_by_key()
    ht = ht.select(spliceai_ds_max=hl.max(ht.values.ds_max))
    ht = ht.annotate_globals(spliceai_version="v1.3")
    ht = ht.checkpoint(output_path, overwrite=True)
    return ht


def run(cadd_outpath: str, spliceai_outpath: str, tmp_dir: str) -> hl.Table:
    """
    Run the insilico predictors workflow to create CADD and SpliceAI Hail Tables.

    :param cadd_outpath: Output path for the CADD Hail Table.
    :param spliceai_outpath: Output path for the SpliceAI Hail Table.
    :param tmp_dir: Temporary directory for intermediate files.
    :return: Hail Table with CADD and SpliceAI scores.
    """
    cadd_ht = create_cadd_grch38_ht(cadd_outpath, tmp_dir)
    spliceai_ht = create_spliceai_grch38_ht(spliceai_outpath)
