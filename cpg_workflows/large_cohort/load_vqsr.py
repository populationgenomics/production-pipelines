import gzip
import logging

import hail as hl

from cpg_utils.hail_batch import genome_build
from cpg_workflows.utils import can_reuse
from gnomad.utils.annotations import annotate_allele_info
from gnomad.utils.sparse_mt import split_info_annotation


def reheader_vcf(vcf_path: str, out_vcf_header_path: str) -> str:
    logging.info('Reheadering VQSR site-only VCF')
    with gzip.open(vcf_path, 'r') as vcf_file:
        with gzip.open(out_vcf_header_path, 'w') as out_file:
            for line in vcf_file:
                if not line.startswith(b'#'):
                    break
                if line.startswith(b'##INFO=<ID=SB'):
                    # rewrite the line with the correct type
                    line = line.replace(b'Number=1', b'Number=.')
                    line = line.replace(b'Type=Float', b'Type=Integer')
                out_file.write(line)
    return out_vcf_header_path


def run(
    pre_vcf_adjusted_ht_path: str,
    vqsr_siteonly_vcf_path: str,
    out_ht_path: str,
):
    load_vqsr(pre_vcf_adjusted_ht_path, vqsr_siteonly_vcf_path, out_ht_path)


def split_info(info_ht: hl.Table) -> hl.Table:
    """
    Splits multi-allelic sites in the provided info Table.

    This function is adapted from `gnomad_methods` (`gnomad.utils.sparse_mt`) to handle
    multi-allelic sites in the info Table. The `AS_lowqual` annotation splitting is
    omitted as it is not used in the VQSR pipeline.

    The `annotate_allele_info` function from `gnomad_methods` is used to split
    multi-allelic sites before splitting the `info` annotation. This ensures that
    all sites in the returned Table are properly annotated with allele-specific
    information.

    Parameters
    ----------
    info_ht : hl.Table
        Info Table containing unsplit multi-allelic sites.

    Returns
    -------
    hl.Table
        Info Table with multi-allelic sites split into bi-allelic sites.
    """
    info_ht = annotate_allele_info(info_ht)
    info_annotations_to_split = ['info']
    info_ht = info_ht.annotate(
        **{
            a: info_ht[a].annotate(
                **split_info_annotation(info_ht[a], info_ht.a_index),
            )
            for a in info_annotations_to_split
        },
        # Vqsr info_ht has no field AS_lowqual
        # AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )

    return info_ht


def load_vqsr(
    pre_vcf_adjusted_ht_path: str,
    vqsr_siteonly_vcf_path: str,
    out_ht_path: str | None = None,
) -> hl.Table:
    """
    Convert VQSR VCF to HT
    """
    if can_reuse(out_ht_path):
        return hl.read_table(str(out_ht_path))

    pre_vcf_adjusted_ht = hl.read_table(str(pre_vcf_adjusted_ht_path))

    reheadered_header_path: str = reheader_vcf(
        vqsr_siteonly_vcf_path,
        'tmp.vcf.gz',
    )

    logging.info(f'AS-VQSR: importing annotations from a site-only VCF {vqsr_siteonly_vcf_path}')
    ht_unsplit = hl.import_vcf(
        str(vqsr_siteonly_vcf_path),
        reference_genome=genome_build(),
        force_bgz=True,
        header_file=reheadered_header_path,
        array_elements_required=False,
    ).rows()

    # VCF has SB fields as float in header:
    # > ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
    # Even though they are lists of ints, e.g. SB=6,11,2,0
    # Hail would fail to parse it, throwing:
    # > java.lang.NumberFormatException: For input string: "6,11,2,0"
    # To mitigate this, we can drop the SB field before the HT is (lazily) parsed.
    # In order words, dropping it before calling ht.write() makes sure that Hail would
    # never attempt to actually parse it.
    # if 'SB' in ht_unsplit.info:
    #     ht_unsplit = ht_unsplit.annotate(info=ht_unsplit.info.drop('SB'))

    # Replace AS_SB_TABLE field in vqsr vcf with correctly formatted array<array<int32>> dtype
    # This field then gets flattened during splitting, eg:
    # eg. chr1:10329	["AC","A"]  [[3,11],[5,19]]	--> [3,11,5,19]
    ht_unsplit = ht_unsplit.annotate(
        info=ht_unsplit.info.annotate(
            AS_SB_TABLE=pre_vcf_adjusted_ht[ht_unsplit.key].info.AS_SB_TABLE,
        ),
    )

    unsplit_count = ht_unsplit.count()

    ht_split = split_info(ht_unsplit)

    # AS_VQSLOD is a string in the VCF, and needs to be converted to float for browser loading
    ht_split = ht_split.annotate(
        info=ht_split.info.annotate(
            AS_VQSLOD=hl.float64(ht_split.info.AS_VQSLOD),
        ),
    )

    split_count = ht_split.count()

    if out_ht_path:
        ht_split.write(str(out_ht_path), overwrite=True)
        ht_split = hl.read_table(str(out_ht_path))
        logging.info(f'Wrote split HT to {out_ht_path}')
    logging.info(f'Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations')
    return ht_split
