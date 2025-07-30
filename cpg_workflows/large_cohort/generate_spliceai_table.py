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
    snvs_path = reference_path('vep_110_mount') + '/spliceai_scores.raw.snv.hg38.vcf.gz'
    indels_path = reference_path('vep_110_mount') + '/spliceai_scores.raw.indel.hg38.vcf.gz'
    # v3_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v3_indel.hg38.vcf.bgz"
    # v4_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v4_new_indels.hg38.vcf.bgz"
    # v3_v4_unscored_path = (
    #     "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v3_v4_unscored_indels.hg38.vcf.bgz"
    # )
    # header_file_path = "gs://gnomad-insilico/spliceai/spliceai.vcf.header"

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
    # v3_indels = import_spliceai_vcf(v3_indels_path)
    # v4_indels = import_spliceai_vcf(v4_indels_path)
    # v3_v4_unscored_indels = import_spliceai_vcf(v3_v4_unscored_path)

    ht = spliceai_snvs.union(
        spliceai_indels,
        # v3_indels,
        # v4_indels,
        # v3_v4_unscored_indels,
    )

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
