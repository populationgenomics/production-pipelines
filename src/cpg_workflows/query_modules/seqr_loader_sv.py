"""
Hail Query functions for seqr loader; SV edition.
"""

import gzip
import logging

import requests

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config, reference_path
from cpg_utils.hail_batch import genome_build
from cpg_workflows.utils import checkpoint_hail, read_hail

# I'm just going to go ahead and steal these constants from their seqr loader
BOTHSIDES_SUPPORT = 'BOTHSIDES_SUPPORT'
GENE_SYMBOL = 'gene_symbol'
GENE_ID = 'gene_id'
MAJOR_CONSEQUENCE = 'major_consequence'
PASS = 'PASS'

# Used to filter mt.info fields.
CONSEQ_PREDICTED_PREFIX = 'PREDICTED_'
NON_GENE_PREDICTIONS = {
    'PREDICTED_INTERGENIC',
    'PREDICTED_NONCODING_BREAKPOINT',
    'PREDICTED_NONCODING_SPAN',
}

# path for downloading this file
GENCODE_GTF_URL = (
    'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/'
    'release_{gencode_release}/gencode.v{gencode_release}.annotation.gtf.gz'
)

PREVIOUS_GENOTYPE_N_ALT_ALLELES = hl.dict(
    {
        # Map of concordance string -> previous n_alt_alleles()
        # Concordant
        frozenset(['TN']): 0,  # 0/0 -> 0/0
        frozenset(['TP']): 2,  # 1/1 -> 1/1
        frozenset(['TN', 'TP']): 1,  # 0/1 -> 0/1
        # Novel
        frozenset(['FP']): 0,  # 0/0 -> 1/1
        frozenset(['TN', 'FP']): 0,  # 0/0 -> 0/1
        # Absent
        frozenset(['FN']): 2,  # 1/1 -> 0/0
        frozenset(['TN', 'FN']): 1,  # 0/1 -> 0/0
        # Discordant
        frozenset(['FP', 'TP']): 1,  # 0/1 -> 1/1
        frozenset(['FN', 'TP']): 2,  # 1/1 -> 0/1
    },
)

GENCODE_FILE_HEADER = [
    'chrom',
    'source',
    'feature_type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'info',
]


# yoinking some methods out of hail_scripts.computed_fields
# removes dependency on submodule completely
def get_expr_for_contig_number(locus: hl.LocusExpression) -> hl.Int32Expression:
    """Convert contig name to contig number"""
    return hl.bind(
        lambda contig: (
            hl.case().when(contig == 'X', 23).when(contig == 'Y', 24).when(contig[0] == 'M', 25).default(hl.int(contig))
        ),
        locus.contig.replace('^chr', ''),
    )


def get_expr_for_xpos(
    locus: hl.LocusExpression | hl.StructExpression,
) -> hl.Int64Expression:
    """Genomic position represented as a single number = contig_number * 10**9 + position.
    This represents chrom:pos more compactly and allows for easier sorting.
    """
    contig_number = get_expr_for_contig_number(locus)
    return hl.int64(contig_number) * 1_000_000_000 + locus.position


def get_cpx_interval(x):
    """
    an example format of CPX_INTERVALS is "DUP_chr1:1499897-1499974"
    Args:
        x (hl.StringExpression instance):
    Returns:
        Struct of CPX components
    """
    type_chr = x.split('_chr')
    chr_pos = type_chr[1].split(':')
    pos = chr_pos[1].split('-')
    return hl.struct(type=type_chr[0], chrom=chr_pos[0], start=hl.int32(pos[0]), end=hl.int32(pos[1]))


def download_gencode_gene_id_mapping(gencode_release: str) -> str:
    """
    This is an inefficient stand-in for now. Swap this out with a more permanent
    location for this file in the resources bucket
    Not suuuuuuper keen on storing a pickled representation, but a minimised JSON
    would be a good middle ground. Parsed ~= 62000 {str: str}
    Args:
        gencode_release (str | int): Which gencode release do you want?
    Returns:
        str - path to localised GTF file
    """

    # this is the thing
    local_gtf = 'localfile.gtf.gz'

    # check for presence in config
    if (config_gtf := get_config().get('gencode_gtf_path')) is not None:
        # treat as GCP
        if config_gtf.startswith('gs://'):
            to_path(config_gtf).copy(local_gtf)
            gtf_path = local_gtf
        else:
            gtf_path = config_gtf
    else:
        gtf_path = GENCODE_GTF_URL.format(gencode_release=gencode_release)

    # if it wasn't a local file, download it
    if gtf_path.startswith(('http://', 'https://')):
        gz_stream = requests.get(gtf_path, stream=True)
        with open(local_gtf, 'wb') as f:
            f.writelines(gz_stream)
        gz_stream.close()

    # OK - I think that covers all the viable situations for now
    return local_gtf


def parse_gtf_from_local(gtf_path: str) -> hl.dict:
    """
    Read over the localised file and read into a dict
    Args:
        gtf_path ():
    Returns:
        the gene lookup dictionary as a Hail DictExpression
    """
    gene_id_mapping = {}
    logging.info(f'Loading {gtf_path}')
    with gzip.open(gtf_path, 'rt') as gencode_file:
        # iterate over this file and do all the things
        for i, line in enumerate(gencode_file):
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != len(GENCODE_FILE_HEADER):
                raise ValueError(f'Unexpected number of fields on line #{i}: {fields}')
            record = dict(zip(GENCODE_FILE_HEADER, fields))
            if record['feature_type'] != 'gene':
                continue
            # parse info field
            info_fields_list = [x.strip().split() for x in record['info'].split(';') if x != '']
            info_fields = {k: v.strip('"') for k, v in info_fields_list}
            gene_id_mapping[info_fields['gene_name']] = info_fields['gene_id'].split('.')[0]
    logging.info('Completed ingestion of gene-ID mapping')
    return hl.literal(gene_id_mapping)


def annotate_cohort_sv(vcf_path: str, out_mt_path: str, checkpoint_prefix: str | None = None):
    """
    Translate an annotated SV VCF into a Seqr-ready format
    Relevant gCNV specific schema
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py
    Relevant gCNV loader script
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/seqr_sv_loading.py
    Args:
        vcf_path (str): Where is the VCF??
        out_mt_path (str): And where do you need output!?
        checkpoint_prefix (str): CHECKPOINT!@!!
    """

    logger = logging.getLogger('annotate_cohort_sv')
    logger.setLevel(logging.INFO)

    logger.info(f'Importing SV VCF {vcf_path}')
    mt = hl.import_vcf(
        vcf_path,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
    )

    # add attributes required for Seqr
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        datasetType='SV',
    )
    if sequencing_type := get_config()['workflow'].get('sequencing_type'):
        # Map to Seqr-style string
        # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        mt = mt.annotate_globals(
            sampleType={
                'genome': 'WGS',
                'exome': 'WES',
                'single_cell': 'RNA',
            }.get(sequencing_type, ''),
        )

    # reimplementation of
    # github.com/populationgenomics/seqr-loading-pipelines..luigi_pipeline/lib/model/sv_mt_schema.py
    mt = mt.annotate_rows(
        sc=mt.info.AC[0],
        sf=mt.info.AF[0],
        sn=mt.info.AN,
        end=mt.info.END,
        end_locus=hl.if_else(
            hl.is_defined(mt.info.END2),
            hl.struct(contig=mt.info.CHR2, position=mt.info.END2),
            hl.struct(contig=mt.locus.contig, position=mt.info.END),
        ),
        sv_callset_Het=mt.info.N_HET,
        sv_callset_Hom=mt.info.N_HOMALT,
        gnomad_svs_ID=mt.info['gnomad_v2.1_sv_SVID'],
        gnomad_svs_AF=mt.info['gnomad_v2.1_sv_AF'],
        gnomad_svs_AC=hl.missing('float64'),
        gnomad_svs_AN=hl.missing('float64'),
        # I DON'T HAVE THESE ANNOTATIONS
        # gnomad_svs_AC=unsafe_cast_int32(mt.info.gnomAD_V2_AC),
        # gnomad_svs_AN=unsafe_cast_int32(mt.info.gnomAD_V2_AN),
        StrVCTVRE_score=hl.parse_float(mt.info.StrVCTVRE),
        filters=hl.or_missing(  # hopefully this plays nicely
            (mt.filters.filter(lambda x: (x != PASS) & (x != BOTHSIDES_SUPPORT))).size() > 0,
            mt.filters,
        ),
        bothsides_support=mt.filters.any(lambda x: x == BOTHSIDES_SUPPORT),
        algorithms=mt.info.ALGORITHMS,
        cpx_intervals=hl.or_missing(
            hl.is_defined(mt.info.CPX_INTERVALS),
            mt.info.CPX_INTERVALS.map(lambda x: get_cpx_interval(x)),
        ),
        sv_types=mt.alleles[1].replace('[<>]', '').split(':', 2),
    )

    # save those changes
    mt = checkpoint_hail(mt, 'initial_annotation_round.mt', checkpoint_prefix)

    # get the Gene-Symbol mapping dict
    gene_id_mapping_file = download_gencode_gene_id_mapping(get_config().get('gencode_release', '42'))
    gene_id_mapping = parse_gtf_from_local(gene_id_mapping_file)

    # OK, NOW IT'S BUSINESS TIME
    conseq_predicted_gene_cols = [
        gene_col
        for gene_col in mt.info
        if (gene_col.startswith(CONSEQ_PREDICTED_PREFIX) and gene_col not in NON_GENE_PREDICTIONS)
    ]

    # register a chain file
    liftover_path = reference_path('liftover_38_to_37')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(str(liftover_path), rg37)

    # annotate with mapped genes
    # Note I'm adding a Flake8 noqa for B023 (loop variable gene_col unbound)
    # I've experimented in a notebook and this seems to perform as expected
    # The homologous small variant seqr_loader method performs a similar function
    # but in a slightly more complicated way (mediated by a method in the S-L-P
    # codebase, so as not to trigger Flake8 evaluation)
    # pos/contig/xpos sourced from
    # seqr-loading-pipelines...luigi_pipeline/lib/model/seqr_mt_schema.py#L12
    mt = mt.annotate_rows(
        sortedTranscriptConsequences=hl.filter(
            hl.is_defined,
            [
                mt.info[gene_col].map(
                    lambda gene: hl.struct(
                        **{
                            GENE_SYMBOL: gene,
                            GENE_ID: gene_id_mapping.get(gene, hl.missing(hl.tstr)),
                            MAJOR_CONSEQUENCE: gene_col.replace(CONSEQ_PREDICTED_PREFIX, '', 1),  # noqa: B023
                        },
                    ),
                )
                for gene_col in conseq_predicted_gene_cols
            ],
        ).flatmap(lambda x: x),
        contig=mt.locus.contig.replace('^chr', ''),
        start=mt.locus.position,
        pos=mt.locus.position,
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(mt.end_locus),
        rg37_locus=hl.liftover(mt.locus, 'GRCh37'),
        rg37_locus_end=hl.or_missing(
            mt.end_locus.position <= hl.literal(hl.get_reference('GRCh38').lengths)[mt.end_locus.contig],
            hl.liftover(
                hl.locus(
                    mt.end_locus.contig,
                    mt.end_locus.position,
                    reference_genome='GRCh38',
                ),
                'GRCh37',
            ),
        ),
        svType=mt.sv_types[0],
        sv_type_detail=hl.if_else(
            mt.sv_types[0] == 'CPX',
            mt.info.CPX_TYPE,
            hl.or_missing(
                (mt.sv_types[0] == 'INS') & (hl.len(mt.sv_types) > 1),
                mt.sv_types[1],
            ),
        ),
        variantId=mt.rsid,
        docId=mt.rsid[0:512],
    )

    # and some more annotation stuff
    mt = mt.annotate_rows(
        transcriptConsequenceTerms=hl.set(
            mt.sortedTranscriptConsequences.map(lambda x: x[MAJOR_CONSEQUENCE]).extend([mt.sv_types[0]]),
        ),
        geneIds=hl.set(
            mt.sortedTranscriptConsequences.filter(lambda x: x[MAJOR_CONSEQUENCE] != 'NEAREST_TSS').map(
                lambda x: x[GENE_ID],
            ),
        ),
        rsid=hl.missing('tstr'),
    )

    # write this output
    mt.write(out_mt_path, overwrite=True)


def annotate_dataset_sv(mt_path: str, out_mt_path: str):
    """
    load the stuff specific to samples in this dataset
    do this after subsetting to specific samples

    Removing the current logic around comparing genotypes to a previous
    callset - doesn't fit with the current implementation

    Args:
        mt_path (str): path to the annotated MatrixTable
        out_mt_path (str): and where do you want it to end up?
    """

    logging.info('Annotating genotypes')

    mt = read_hail(mt_path)
    is_called = hl.is_defined(mt.GT)
    num_alt = hl.if_else(is_called, mt.GT.n_alt_alleles(), -1)
    # was_previously_called = hl.is_defined(mt.CONC_ST) & ~mt.CONC_ST.contains('EMPTY')
    # prev_num_alt = hl.if_else(
    #     was_previously_called, PREVIOUS_GENOTYPE_N_ALT_ALLELES[hl.set(mt.CONC_ST)], -1
    # )
    # concordant_genotype = num_alt == prev_num_alt
    # discordant_genotype = (num_alt != prev_num_alt) & (prev_num_alt > 0)
    # novel_genotype = (num_alt != prev_num_alt) & (prev_num_alt == 0)
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.GQ,
                cn=mt.RD_CN,
                num_alt=num_alt,
                # prev_num_alt=hl.or_missing(discordant_genotype, prev_num_alt),
                # prev_call=hl.or_missing(
                #     is_called, was_previously_called & concordant_genotype
                # ),
                # new_call=hl.or_missing(
                #     is_called, ~was_previously_called | novel_genotype
                # ),
            ),
        ),
    )

    def _genotype_filter_samples(fn):
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

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
    def _filter_num_alt(i, g):
        return i == g.num_alt

    @_capture_i_decorator
    def _filter_samples_gq(i, g):
        return (g.gq >= i) & (g.gq < i + 10)

    @_capture_i_decorator
    def _filter_sample_cn(i, g):
        return g.cn == i

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/sv_mt_schema.py#L221
    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/seqr_mt_schema.py#L251
    mt = mt.annotate_rows(
        # omit samples field for GATKSV callsets. Leaving this here as likely needed
        # for gCNV specific callsets (maybe)
        # samples=_genotype_filter_samples(lambda g: True),
        # samples_new_alt=_genotype_filter_samples(
        #     lambda g: g.new_call | hl.is_defined(g.prev_num_alt)
        # ),
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{'%i' % i: _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq_sv=hl.struct(
            **{('%i_to_%i' % (i, i + 10)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 90, 10)},
        ),
        # As per `samples` field, I beleive CN stats should only be generated for gCNV only
        # callsets. In particular samples_cn_2 is used to select ALT_ALT variants,
        # presumably because this genotype is only asigned when this CN is alt (ie on chrX)
        # samples_cn=hl.struct(
        #     **{
        #         f'{i}': _genotype_filter_samples(_filter_sample_cn(i))
        #         for i in range(0, 4, 1)
        #     },
        #     gte_4=_genotype_filter_samples(lambda g: g.cn >= 4),
        # ),
    )

    logging.info('Genotype fields annotated')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Written  SV MT to {out_mt_path}')
