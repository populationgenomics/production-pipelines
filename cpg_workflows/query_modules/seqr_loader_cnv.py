"""
Hail Query functions for seqr loader; SV edition.
"""
import datetime
import logging

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build, reference_path
from cpg_workflows.query_modules.seqr_loader_sv import get_expr_for_xpos, parse_gtf_from_local, download_gencode_gene_id_mapping
from cpg_workflows.utils import read_hail, checkpoint_hail


# I'm just going to go ahead and steal these constants from their seqr loader
GENE_SYMBOL = 'gene_symbol'
GENE_ID = 'gene_id'
MAJOR_CONSEQUENCE = 'major_consequence'

# Used to filter mt.info fields.
CONSEQ_PREDICTED_PREFIX = 'PREDICTED_'
NON_GENE_PREDICTIONS = {
    'PREDICTED_INTERGENIC',
    'PREDICTED_NONCODING_BREAKPOINT',
    'PREDICTED_NONCODING_SPAN',
}


def parse_genes(gene_col: hl.expr.StringExpression) -> hl.expr.SetExpression:
    """
    Convert a string-ified gene list to a set()
    """
    return hl.set(gene_col.split(',').filter(
        lambda gene: ~hl.set({'None', 'null', 'NA', ''}).contains(gene)
    ).map(
        lambda gene: gene.split(r'\.')[0]
    ))


def hl_agg_collect_set_union(gene_col: hl.expr.SetExpression) -> hl.expr.SetExpression:
    return hl.flatten(hl.agg.collect_as_set(gene_col))

def annotate_cohort_gcnv(
    vcf_path: str, out_mt_path: str, checkpoint_prefix: str | None = None
):
    """
    Translate an annotated gCNV VCF into a Seqr-ready format
    Relevant gCNV specific schema
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/gcnv_mt_schema.py
    Relevant gCNV loader script
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/seqr_gcnv_loading.py
    Args:
        vcf_path (str): Where is the VCF??
        out_mt_path (str): And where do you need output!?
        checkpoint_prefix (str): CHECKPOINT!@!!
    """

    logger = logging.getLogger('annotate_cohort_gcnv')
    logger.setLevel(logging.INFO)

    logger.info(f'Importing SV VCF {vcf_path}')
    mt = hl.import_vcf(
        vcf_path,
        array_elements_required=False,
        force_bgz=True,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
    )

    # add attributes required for Seqr
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        datasetType='CNV',
        sampleType='WES'
    )

    # apply variant_qc annotations
    mt = hl.variant_qc(mt)

    # reimplementation of
    # github.com/populationgenomics/seqr-loading-pipelines..luigi_pipeline/lib/model/sv_mt_schema.py
    mt = mt.annotate_rows(
        sc=mt.variant_qc.AC[0],
        sf=mt.variant_qc.AF[0],
        sn=mt.variant_qc.AN,
        end=mt.info.END,
        sv_callset_Het=mt.info.N_HET,
        sv_callset_Hom=mt.info.N_HOMALT,
        gnomad_svs_ID=mt.info['gnomad_v2.1_sv_SVID'],
        gnomad_svs_AF=mt.info['gnomad_v2.1_sv_AF'],
        gnomad_svs_AC=hl.missing('float64'),
        gnomad_svs_AN=hl.missing('float64'),
        StrVCTVRE_score=hl.parse_float(mt.info.StrVCTVRE),
        svType=mt.alleles[1].replace('[<>]', ''),
        start=mt.locus.position,
        pos=mt.locus.position,
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(hl.struct(contig=mt.locus.contig, position=mt.info.END)),
    )

    # save those changes
    mt = checkpoint_hail(mt, 'initial_annotation_round.mt', checkpoint_prefix)

    # get the Gene-Symbol mapping dict
    gene_id_mapping_file = download_gencode_gene_id_mapping(
        get_config().get('gencode_release', '42')
    )
    gene_id_mapping = parse_gtf_from_local(gene_id_mapping_file)

    # find all the column names which contain Gene symbols
    conseq_predicted_gene_cols = [
        gene_col
        for gene_col in mt.info
        if (
                gene_col.startswith(CONSEQ_PREDICTED_PREFIX)
                and gene_col not in NON_GENE_PREDICTIONS
        )
    ]

    # bank all those symbols before overwriting them - may not be required
    # might have to drop this for the Seqr ingest
    mt = mt.annotate_rows(
        geneSymbols=hl.set(hl.filter(
            hl.is_defined,
            [
                mt.info[gene_col] for gene_col in conseq_predicted_gene_cols
            ]
        ).flatmap(lambda x: x))
    )

    # overwrite symbols with ENSG IDs in these columns
    # not sure why this is required, I think SV annotation came out
    # with ENSGs from the jump, but this is all symbols
    for col_name in conseq_predicted_gene_cols:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **{
                    col_name: hl.map(lambda gene: gene_id_mapping.get(gene, gene), mt.info[col_name])
                }
            )
        )

    mt = mt.annotate_rows(
        # this expected mt.variant_name to be present, and it's not
        variantId=hl.format(f"%s_%s_{datetime.date.today():%m%d%Y}", mt.rsid, mt.svType),
        geneIds=hl.set(hl.filter(
            hl.is_defined,
            [
                mt.info[gene_col] for gene_col in conseq_predicted_gene_cols
            ]
        ).flatmap(lambda x: x))
    )

    lof_genes = hl.set(mt.info.PREDICTED_LOF)
    major_consequence_genes = lof_genes | hl.set(mt.info.PREDICTED_COPY_GAIN)
    mt = mt.annotate_rows(
        sortedTranscriptConsequences=hl.map(
            lambda gene: hl.Struct(
                gene_id=gene,
                major_consequence=hl.or_missing(
                    major_consequence_genes.contains(gene),
                    hl.if_else(
                        lof_genes.contains(gene),
                        "LOF",
                        "COPY_GAIN",
                    ),
                ),
            ),
            mt.geneIds,
        )
    )

    # transcriptConsequenceTerms
    default_consequences = [hl.format('gCNV_%s', mt.svType)]
    gene_major_consequences = hl.array(hl.set(
        mt.sortedTranscriptConsequences
        .filter(lambda x: hl.is_defined(x.major_consequence))
        .map(lambda x: x.major_consequence)
    ))
    mt = mt.annotate_rows(
        transcriptConsequenceTerms=gene_major_consequences.extend(default_consequences),
        docId=mt.variantId[0:512],
    )

    # write this output
    mt.write(out_mt_path, overwrite=True)
