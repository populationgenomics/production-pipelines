import hail as hl
import logging

from cpg_utils import Path
from cpg_utils.hail_batch import genome_build
from cpg_workflows.utils import can_reuse


def run(
    site_only_vcf_path: Path,
    out_ht_path: Path,
):
    load_vqsr(site_only_vcf_path, out_ht_path)


def load_vqsr(
    site_only_vcf_path: Path,
    out_ht_path: Path | None = None,
) -> hl.Table:
    """
    Convert VQSR VCF to HT
    """
    if can_reuse(out_ht_path):
        return hl.read_table(str(out_ht_path))

    logging.info(
        f'AS-VQSR: importing annotations from a site-only VCF {site_only_vcf_path}'
    )
    ht = hl.import_vcf(
        str(site_only_vcf_path),
        reference_genome=genome_build(),
        force_bgz=True,
    ).rows()

    # VCF has SB fields as float in header:
    # > ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
    # Even though they are lists of ints, e.g. SB=6,11,2,0
    # Hail would fail to parse it, throwing:
    # > java.lang.NumberFormatException: For input string: "6,11,2,0"
    # To mitigate this, we can drop the SB field before the HT is (lazily) parsed.
    # In order words, dropping it before calling ht.write() makes sure that Hail would
    # never attempt to actually parse it.
    if 'SB' in ht.info:
        ht = ht.annotate(info=ht.info.drop('SB'))

    # Dropping also all INFO/AS* annotations as well as InbreedingCoeff, as they are
    # causing problems splitting multiallelics after parsing by Hail, when Hail attempts
    # to subset them by allele index. For example, for these 2 variants:
    # chr1    10145   .       AAC     A,TAC   .       VQSRTrancheINDEL99.50to99.90    AC=0,0;AC_raw=1,1;AS_FS=0,0;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=45.5636,46.0067;AS_MQRankSum=0.092,1.34;AS_QD=3.64286,1;AS_QUALapprox=51,20;AS_ReadPosRankSum=0.657,1.128;AS_SB_TABLE=15,15|1,1|1,1;AS_SOR=0.693147,0.693147;AS_VQSLOD=-1.9389;AS_VarDP=14,20;AS_culprit=AS_MQRankSum;DP=1908;FS=0;MQ=45.7857;MQRankSum=1.34;NEGATIVE_TRAIN_SITE;QD=2.08824;QUALapprox=71;ReadPosRankSum=1.128;SB=15,15,2,2;SOR=0.693147;VarDP=34
    # chr1    10146   .       AC      A       .       VQSRTrancheINDEL99.50to99.90    AC=4;AC_raw=5;AS_FS=5.75068;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=41.3793;AS_MQRankSum=1.033;AS_QD=12.1209;AS_QUALapprox=1103;AS_ReadPosRankSum=-0.875;AS_SB_TABLE=18,12|28,33;AS_SOR=0.611231;AS_VQSLOD=-2.0660;AS_VarDP=91;AS_culprit=AS_MQRankSum;DP=1727;FS=5.75068;MQ=41.3793;MQRankSum=1.033;NEGATIVE_TRAIN_SITE;QD=12.1209;QUALapprox=1103;ReadPosRankSum=-0.875;SB=18,12,28,33;SOR=0.611231;VarDP=91
    # The first one has 2 alleles, and one AS_FilterStatus value, same as the second
    # one, with one AS_FilterStatus value. So for the first one indexing by allele
    # index would work, but for the second one it would throw an index out of bounds:
    # `HailException: array index out of bounds: index=1, length=1`
    ht = ht.annotate(info=ht.info.drop(*[f for f in ht.info if f.startswith('AS_')]))

    unsplit_count = ht.count()
    ht = hl.split_multi_hts(ht)
    if out_ht_path:
        ht.write(str(out_ht_path), overwrite=True)
        ht = hl.read_table(str(out_ht_path))
        logging.info(f'Wrote split HT to {out_ht_path}')
    split_count = ht.count()
    logging.info(
        f'Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations'
    )
    return ht
