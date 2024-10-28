"""
Hail Query functions for seqr loader; CNV edition.
"""

import datetime
import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils.hail_batch import genome_build, init_batch
from cpg_workflows.query_modules.seqr_loader_sv import get_expr_for_xpos, parse_gtf_from_local

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


def annotate_cohort_gcnv(vcf: str, mt_out: str, gencode: str, checkpoint: str, *args, **kwargs):
    """
    Translate an annotated gCNV VCF into a Seqr-ready format
    Relevant gCNV specific schema
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/gcnv_mt_schema.py
    Relevant gCNV loader script
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/seqr_gcnv_loading.py
    Args:
        vcf (str): Where is the VCF??
        mt_out (str): And where do you need output!?
        gencode (str): The path to a compressed GENCODE GTF file
        checkpoint (str): location we can write checkpoints to
    """

    logging.info(f'Importing SV VCF {vcf}')
    mt = hl.import_vcf(
        vcf,
        array_elements_required=False,
        force_bgz=True,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
    )

    # add attributes required for Seqr
    mt = mt.annotate_globals(
        sourceFilePath=vcf,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        datasetType='SV',
        sampleType='WES',
    )

    # reimplementation of
    # github.com/populationgenomics/seqr-loading-pipelines..luigi_pipeline/lib/model/gcnv_mt_schema.py
    mt = mt.annotate_rows(
        contig=mt.locus.contig.replace('^chr', ''),
        start=mt.locus.position,
        pos=mt.locus.position,
        # todo @MattWellie - review use of AC_Orig vs. AC (post-qc)
        sc=mt.info.AC_Orig[0],
        sf=mt.info.AF_Orig[0],
        sn=mt.info.AN_Orig,
        end=mt.info.END,
        sv_callset_Het=mt.info.N_HET,
        sv_callset_Hom=mt.info.N_HOMALT,
        gnomad_svs_ID=mt.info['gnomad_v2.1_sv_SVID'],
        gnomad_svs_AF=mt.info['gnomad_v2.1_sv_AF'],
        gnomad_svs_AC=hl.missing('float64'),
        gnomad_svs_AN=hl.missing('float64'),
        StrVCTVRE_score=hl.parse_float(mt.info.StrVCTVRE),
        svType=mt.alleles[1].replace('[<>]', ''),
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(hl.struct(contig=mt.locus.contig, position=mt.info.END)),
        num_exon=hl.agg.max(mt.NP),
    )

    # find all the column names which contain Gene symbols
    conseq_predicted_gene_cols = [
        gene_col
        for gene_col in mt.info
        if (gene_col.startswith(CONSEQ_PREDICTED_PREFIX) and gene_col not in NON_GENE_PREDICTIONS)
    ]

    # bank all those symbols before overwriting them - may not be required
    # might have to drop this for the Seqr ingest
    mt = mt.annotate_rows(
        geneSymbols=hl.set(
            hl.filter(
                hl.is_defined,
                [mt.info[gene_col] for gene_col in conseq_predicted_gene_cols],
            ).flatmap(
                lambda x: x,
            ),
        ),
    )

    mt = mt.checkpoint(join(checkpoint, 'pre-gene_annotation.mt'))

    # this next section is currently failing - the dictionary of genes is too large
    # to be used in an annotation expression. At least... I think it is
    # for i, chunks in enumerate(chunks(gene_id_mapping, 100)):

    # get the Gene-Symbol mapping dict
    gene_id_mappings = parse_gtf_from_local(gencode, chunk_size=15000)

    # overwrite symbols with ENSG IDs in these columns
    # not sure why this is required, I think SV annotation came out
    # with ENSGs from the jump, but this is all symbols
    for i, gene_id_mapping in enumerate(gene_id_mappings):
        logging.info(f'Processing gene ID mapping chunk {i}')
        for col_name in conseq_predicted_gene_cols:
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **{col_name: hl.map(lambda gene: gene_id_mapping.get(gene, gene), mt.info[col_name])},
                ),
            )

        mt = mt.checkpoint(join(checkpoint, f'fragment_{i}.mt'))

    mt = mt.annotate_rows(
        # this expected mt.variant_name to be present, and it's not
        variantId=hl.format(f'%s_%s_{datetime.date.today():%m%d%Y}', mt.rsid, mt.svType),
        geneIds=hl.set(
            hl.filter(
                hl.is_defined,
                [mt.info[gene_col] for gene_col in conseq_predicted_gene_cols],
            ).flatmap(
                lambda x: x,
            ),
        ),
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
                        'LOF',
                        'COPY_GAIN',
                    ),
                ),
            ),
            mt.geneIds,
        ),
    )

    # transcriptConsequenceTerms
    default_consequences = [hl.format('gCNV_%s', mt.svType)]
    gene_major_consequences = hl.array(
        hl.set(
            mt.sortedTranscriptConsequences.filter(lambda x: hl.is_defined(x.major_consequence)).map(
                lambda x: x.major_consequence,
            ),
        ),
    )
    mt = mt.annotate_rows(
        transcriptConsequenceTerms=gene_major_consequences.extend(default_consequences),
        docId=mt.variantId[0:512],
    )

    # write this output
    mt.write(mt_out, overwrite=True)


def annotate_dataset_gcnv(mt_in: str, mt_out: str, *args, **kwargs):
    """
    process data specific to samples in this dataset
    do this after sub-setting to specific samples
    Args:
        mt_in (str): path to the annotated MatrixTable
        mt_out (str): and where do you want it to end up?
    """

    logging.info('Annotating genotypes')

    mt = hl.read_matrix_table(mt_in)

    # adding in the GT here, that may cause problems later?
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.QA,
                cn=mt.CN,
                end=mt.end,
                num_exon=mt.NP,
                start=mt.start,
                geneIds=mt.geneIds,
                gt=mt.GT,
            ),
        ),
    )

    def _genotype_filter_samples(fn):
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    def _genotype_filter_samples_cn2():
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(lambda g: ((g.gt.is_haploid()) & (g.cn == 2))).map(lambda g: g.sample_id))

    # top level - decorator
    def _capture_i_decorator(func):
        # call the returned_function(i) which locks in the value of i
        def _inner_filter(i):
            # the _genotype_filter_samples will call this _func with g
            def _func(g):
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_sample_cn(i, g):
        return g.cn == i

    @_capture_i_decorator
    def _filter_samples_gq(i, g):
        return (g.gq >= i) & (g.gq < i + 10)

    @_capture_i_decorator
    def _filter_num_alt(i, g):
        return i == g.gt.n_alt_alleles()

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/gcnv_mt_schema.py
    mt = mt.annotate_rows(
        # # samples is bad, we do not want it
        # samples=_genotype_filter_samples(lambda g: True),
        # dubious about this annotation - expected field is qs, I'm using gq, derived from CNQ
        samples_qs=hl.struct(
            **{('%i_to_%i' % (i, i + 10)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 1000, 10)},
            **{'gt_1000': _genotype_filter_samples(lambda g: g.gq >= 1000)},
        ),
        # ok, here's what we're
        samples_cn=hl.struct(
            **{str(i): _genotype_filter_samples(_filter_sample_cn(i)) for i in [0, 1, 3]},
            **{'gte_4': _genotype_filter_samples(lambda g: g.cn >= 4)},
            # and a special case for male sex-chrom CN
            **{'2': _genotype_filter_samples_cn2()},
        ),
        # re-embedding the samples_num_alt, derived from hl.call().n_alt_alleles()
        samples_num_alt=hl.struct(**{'%i' % i: _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
    )

    # removing GT again, out of an abundance of caution
    # adding in the GT here, that may cause problems later?
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.QA,
                cn=mt.CN,
                end=mt.end,
                num_exon=mt.NP,
                start=mt.start,
                geneIds=mt.geneIds,
            ),
        ),
    )
    logging.info('Genotype fields annotated')
    mt.describe()
    mt.write(mt_out, overwrite=True)
    logging.info(f'Written gCNV MT to {mt_out}')


def cli_main():
    """
    command line entrypoint
    """
    # enable Info-level logging
    logging.basicConfig(level=logging.INFO)

    init_batch()

    # set up an argument parser to allow two separate entrypoints
    parser = ArgumentParser()
    # these arguments are used for both entrypoints
    parser.add_argument('--mt_out', help='Path to the MatrixTable, input or output', required=True)
    parser.add_argument('--checkpoint', help='Dir to write checkpoints to', required=True)
    subparsers = parser.add_subparsers()

    # a parser for the AnnotateCohort method
    cohort_parser = subparsers.add_parser('cohort')
    cohort_parser.add_argument('--vcf', help='Path to input VCF file')
    cohort_parser.add_argument('--gencode', help='Path to input gencode GTF file')
    cohort_parser.set_defaults(func=annotate_cohort_gcnv)

    # a parser for the AnnotateDataset method
    cohort_parser = subparsers.add_parser('dataset')
    cohort_parser.add_argument('--mt_in', help='Path to input MT')
    cohort_parser.set_defaults(func=annotate_dataset_gcnv)

    args = parser.parse_args()

    args.func(**vars(args))


if __name__ == '__main__':
    cli_main()
