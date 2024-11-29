"""
Standalone script to run the AnnotateCohort process in the seqr_loader workflow.
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils import to_path, Path
from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import genome_build, init_batch

from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import can_reuse, checkpoint_hail, get_logger
from hail_scripts.computed_fields import variant_id, vep


def load_vqsr(vcf_path: str, ht_path: Path) -> hl.Table:
    """
    Convert VQSR VCF to HT, and checkpoints

    Args:
        vcf_path ():
        ht_path ():
    """
    if can_reuse(ht_path):
        get_logger().info(f'Reading VQSR checkpoint from {ht_path}')
        return hl.read_table(str(ht_path))
    get_logger().info(f'AS-VQSR: importing annotations from a site-only VCF {vcf_path}')
    vqsr_ht = hl.import_vcf(
        vcf_path,
        reference_genome=genome_build(),
        force_bgz=True,
    ).rows()

    # two comments copied in from the previous implementation, unsure if these are still valid

    # VCF has SB fields as float in header:
    # > ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
    # Even though they are lists of ints, e.g. SB=6,11,2,0
    # Hail would fail to parse it, throwing:
    # > java.lang.NumberFormatException: For input string: "6,11,2,0"
    # To mitigate this, we can drop the SB field before the HT is (lazily) parsed.
    # In order words, dropping it before calling ht.write() makes sure that Hail would
    # never attempt to actually parse it.

    # Dropping also all INFO/AS* annotations as well as InbreedingCoeff, as they are
    # causing problems splitting multiallelics after parsing by Hail, when Hail attempts
    # to subset them by allele index. For example, for these 2 variants:
    # chr1    10145   .       AAC     A,TAC   .       VQSRTrancheINDEL99.50to99.90    AC=0,0;AC_raw=1,1;AS_FS=0,0;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=45.5636,46.0067;AS_MQRankSum=0.092,1.34;AS_QD=3.64286,1;AS_QUALapprox=51,20;AS_ReadPosRankSum=0.657,1.128;AS_SB_TABLE=15,15|1,1|1,1;AS_SOR=0.693147,0.693147;AS_VQSLOD=-1.9389;AS_VarDP=14,20;AS_culprit=AS_MQRankSum;DP=1908;FS=0;MQ=45.7857;MQRankSum=1.34;NEGATIVE_TRAIN_SITE;QD=2.08824;QUALapprox=71;ReadPosRankSum=1.128;SB=15,15,2,2;SOR=0.693147;VarDP=34
    # chr1    10146   .       AC      A       .       VQSRTrancheINDEL99.50to99.90    AC=4;AC_raw=5;AS_FS=5.75068;AS_FilterStatus=VQSRTrancheINDEL99.50to99.90;AS_MQ=41.3793;AS_MQRankSum=1.033;AS_QD=12.1209;AS_QUALapprox=1103;AS_ReadPosRankSum=-0.875;AS_SB_TABLE=18,12|28,33;AS_SOR=0.611231;AS_VQSLOD=-2.0660;AS_VarDP=91;AS_culprit=AS_MQRankSum;DP=1727;FS=5.75068;MQ=41.3793;MQRankSum=1.033;NEGATIVE_TRAIN_SITE;QD=12.1209;QUALapprox=1103;ReadPosRankSum=-0.875;SB=18,12,28,33;SOR=0.611231;VarDP=91
    # The first one has 2 alleles, and one AS_FilterStatus value, same as the second
    # one, with one AS_FilterStatus value. So for the first one indexing by allele
    # index would work, but for the second one it would throw an index out of bounds:
    # `HailException: array index out of bounds: index=1, length=1`
    vqsr_ht = vqsr_ht.annotate(info=vqsr_ht.info.drop(*[f for f in vqsr_ht.info if (f.startswith('AS_') or f == 'SB')]))
    vqsr_ht = vqsr_ht.checkpoint(str(ht_path), overwrite=True)
    return vqsr_ht


def annotate_cohort(
    mt_path: str,
    out_mt_path: str,
    vep_ht_path: str,
    checkpoint_prefix: str,
    vqsr_vcf_path: str | None = None,
    long_read: bool = False,
) -> None:
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR annotations.

    Args:
        mt_path ():
        out_mt_path ():
        vep_ht_path ():
        checkpoint_prefix ():
        vqsr_vcf_path ():
        long_read ():

    Returns:
        Nothing, but hopefully writes out a new MT
    """

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    # we satisfy this requirement by calling `init_batch()` in the query_command wrapper
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    # hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/permissions.20issues/near/398711114
    # don't override the block size, as it explodes the number of partitions when processing TB+ datasets
    # Each partition comes with some computational overhead, it's to be seen whether the standard block size
    # is viable for QOB in large datasets... Requires a test
    mt = hl.read_matrix_table(mt_path)
    get_logger().info(f'Imported MT from {mt_path} as {mt.n_partitions()} partitions')

    # Annotate VEP. Do it before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(vep_ht_path)
    get_logger().info(
        f'Adding VEP annotations into the Matrix Table from {vep_ht_path}. VEP loaded as {vep_ht.n_partitions()} partitions',
    )
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # in our long-read VCFs, AF is present as an Entry field, so we need to drop it from entries,
    # and then recompute AC/AF/AN correctly from the variant QC table
    # do this prior to splitting multiallelics, as the AF/AC needs to be generated per-original ALT allele
    # currently not an issue, as our long-read VCFs are not multiallelic, but they could be in future
    if long_read:
        mt = mt.drop('AF')
        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                AF=mt.variant_qc.AF,
                AN=mt.variant_qc.AN,
                AC=mt.variant_qc.AC,
            ),
        )
        mt = mt.drop('variant_qc')

    # Splitting multi-allelics. We do not handle AS info fields here - we handle
    # them when loading VQSR instead, and populate entire "info" from VQSR.
    mt = hl.split_multi_hts(mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles))
    mt = checkpoint_hail(mt, 'mt-vep-split.mt', checkpoint_prefix)

    if vqsr_vcf_path:
        get_logger().info('Adding VQSR annotations into the Matrix Table')
        vqsr_checkpoint = to_path(checkpoint_prefix) / 'vqsr.ht'
        vqsr_ht = load_vqsr(vqsr_vcf_path, vqsr_checkpoint)
        mt = mt.annotate_globals(**vqsr_ht.index_globals())
        mt = mt.annotate_rows(
            # vqsr_ht has info annotation split by allele, plus the new AS-VQSR annotations
            info=vqsr_ht[mt.row_key].info,
            filters=mt.filters.union(vqsr_ht[mt.row_key].filters).filter(lambda val: val != 'PASS'),
        )
        mt = checkpoint_hail(mt, 'mt-vep-split-vqsr.mt', checkpoint_prefix)

    ref_ht = hl.read_table(reference_path('seqr_combined_reference_data'))
    clinvar_ht = hl.read_table(reference_path('seqr_clinvar'))

    get_logger().info('Annotating with seqr-loader fields: round 1')

    # split the AC/AF attributes into separate entries, overwriting the array in INFO
    # these elements become a 1-element array
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AF=[mt.info.AF[mt.a_index - 1]],
            AC=[mt.info.AC[mt.a_index - 1]],
        ),
    )

    get_logger().info('Annotating with clinvar and munging annotation fields')
    mt = mt.annotate_rows(
        # still taking just a single value here for downstream compatibility in Seqr
        AC=mt.info.AC[0],
        AF=mt.info.AF[0],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=variant_id.get_expr_for_variant_ids(mt.locus_old, mt.alleles_old),
        sortedTranscriptConsequences=vep.get_expr_for_vep_sorted_transcript_consequences_array(mt.vep),
        variantId=variant_id.get_expr_for_variant_id(mt),
        contig=variant_id.get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=variant_id.get_expr_for_xpos(mt.locus),
        xstart=variant_id.get_expr_for_xpos(mt.locus),
        xstop=variant_id.get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
    )

    # this was previously executed in the MtToEs job, as it wasn't possible on QoB
    get_logger().info('Adding GRCh37 coords')
    liftover_path = reference_path('liftover_38_to_37')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_path, rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))

    # only remove InbreedingCoeff if present (post-VQSR)
    if 'InbreedingCoeff' in mt.info:
        mt = mt.annotate_rows(info=mt.info.drop('InbreedingCoeff'))

    get_logger().info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)',
    )
    mt = mt.annotate_rows(
        domains=vep.get_expr_for_vep_protein_domains_set_from_sorted(mt.sortedTranscriptConsequences),
        transcriptConsequenceTerms=vep.get_expr_for_vep_consequence_terms_set(mt.sortedTranscriptConsequences),
        transcriptIds=vep.get_expr_for_vep_transcript_ids_set(mt.sortedTranscriptConsequences),
        mainTranscript=vep.get_expr_for_worst_transcript_consequence_annotations_struct(
            mt.sortedTranscriptConsequences,
        ),
        geneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences, only_coding_genes=True),
        cadd=mt.ref_data.cadd,
        dbnsfp=mt.ref_data.dbnsfp,
        geno2mp=mt.ref_data.geno2mp,
        gnomad_exomes=mt.ref_data.gnomad_exomes,
        gnomad_exome_coverage=mt.ref_data.gnomad_exome_coverage,
        gnomad_genomes=mt.ref_data.gnomad_genomes,
        gnomad_genome_coverage=mt.ref_data.gnomad_genome_coverage,
        eigen=mt.ref_data.eigen,
        exac=mt.ref_data.exac,
        g1k=mt.ref_data.g1k,
        mpc=mt.ref_data.mpc,
        primate_ai=mt.ref_data.primate_ai,
        splice_ai=mt.ref_data.splice_ai,
        topmed=mt.ref_data.topmed,
        clinvar=hl.struct(
            **{
                'allele_id': mt.clinvar_data.info.ALLELEID,
                'clinical_significance': hl.delimit(mt.clinvar_data.info.CLNSIG),
                'gold_stars': mt.clinvar_data.gold_stars,
            },
        ),
    )
    mt = mt.annotate_globals(
        sourceFilePath=mt_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
    )
    if sequencing_type := config_retrieve(['workflow', 'sequencing_type']):
        # Map to Seqr-style string
        # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        mt = mt.annotate_globals(
            sampleType={
                'genome': 'WGS',
                'exome': 'WES',
                'single_cell': 'RNA',
            }.get(sequencing_type, ''),
        )

    get_logger().info('Done:')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    get_logger().info(f'Written final matrix table into {out_mt_path}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input MatrixTable to annotate')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    parser.add_argument('--vep-ht', required=True, help='HT with VEP annotations')
    parser.add_argument('--checkpoint-prefix', required=True, help='Checkpoint prefix')
    parser.add_argument('--vqsr-vcf', required=False, help='Site-only VQSR VCF')
    parser.add_argument('--long-read', action='store_true', help='Long-read VCF')
    args = parser.parse_args()
    annotate_cohort(
        mt_path=args.input,
        out_mt_path=args.output,
        vep_ht_path=args.vep_ht,
        checkpoint_prefix=args.checkpoint_prefix,
        vqsr_vcf_path=args.site_only_vqsr_vcf,
        long_read=args.long_read,
    )


if __name__ == '__main__':
    cli_main()
