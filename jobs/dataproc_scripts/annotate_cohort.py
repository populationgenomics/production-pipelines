#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import os
import logging

import click
import coloredlogs
import hail as hl

from cpg_utils.hail_batch import reference_path, genome_build
from cpg_utils.config import get_config

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='DEBUG', fmt=fmt)


# Consequence terms in order of severity (more severe to less severe) as estimated by Ensembl.
# See https://ensembl.org/info/genome/variation/prediction/predicted_data.html
CONSEQUENCE_TERMS = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

# hail DictExpression that maps each CONSEQUENCE_TERM to it's rank in the list
CONSEQUENCE_TERM_RANK_LOOKUP = hl.dict(
    {term: rank for rank, term in enumerate(CONSEQUENCE_TERMS)}
)


OMIT_CONSEQUENCE_TERMS = [
    "upstream_gene_variant",
    "downstream_gene_variant",
]


def get_expr_for_vep_consequence_terms_set(vep_transcript_consequences_root):
    return hl.set(
        vep_transcript_consequences_root.flatmap(lambda c: c.consequence_terms)
    )


def get_expr_for_vep_gene_ids_set(
    vep_transcript_consequences_root, only_coding_genes=False
):
    """Expression to compute the set of gene ids in VEP annotations for this variant.
    Args:
        vep_transcript_consequences_root (ArrayExpression): VEP transcript_consequences root in the struct
        only_coding_genes (bool): If set to True, non-coding genes will be excluded.
    Return:
        SetExpression: expression
    """

    expr = vep_transcript_consequences_root

    if only_coding_genes:
        expr = expr.filter(lambda c: hl.or_else(c.biotype, "") == "protein_coding")

    return hl.set(expr.map(lambda c: c.gene_id))


def get_expr_for_vep_protein_domains_set(vep_transcript_consequences_root):
    return hl.set(
        vep_transcript_consequences_root.flatmap(
            lambda c: c.domains.map(lambda domain: domain.db + ":" + domain.name)
        )
    )


PROTEIN_LETTERS_1TO3 = hl.dict(
    {
        "A": "Ala",
        "C": "Cys",
        "D": "Asp",
        "E": "Glu",
        "F": "Phe",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "K": "Lys",
        "L": "Leu",
        "M": "Met",
        "N": "Asn",
        "P": "Pro",
        "Q": "Gln",
        "R": "Arg",
        "S": "Ser",
        "T": "Thr",
        "V": "Val",
        "W": "Trp",
        "Y": "Tyr",
        "X": "Ter",
        "*": "Ter",
        "U": "Sec",
    }
)


HGVSC_CONSEQUENCES = hl.set(
    ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]
)


def get_expr_for_formatted_hgvs(csq):
    return hl.cond(
        hl.is_missing(csq.hgvsp) | HGVSC_CONSEQUENCES.contains(csq.major_consequence),
        csq.hgvsc.split(":")[-1],
        hl.cond(
            csq.hgvsp.contains("=") | csq.hgvsp.contains("%3D"),
            hl.bind(
                lambda protein_letters: "p."
                + protein_letters
                + hl.str(csq.protein_start)
                + protein_letters,
                hl.delimit(
                    csq.amino_acids.split("").map(
                        lambda l: PROTEIN_LETTERS_1TO3.get(l)
                    ),
                    "",
                ),
            ),
            csq.hgvsp.split(":")[-1],
        ),
    )


def get_expr_for_vep_sorted_transcript_consequences_array(
    vep_root, include_coding_annotations=True, omit_consequences=OMIT_CONSEQUENCE_TERMS
):
    """Sort transcripts by 3 properties:
        1. coding > non-coding
        2. transcript consequence severity
        3. canonical > non-canonical
    so that the 1st array entry will be for the coding, most-severe, canonical transcript (assuming
    one exists).
    Also, for each transcript in the array, computes these additional fields:
        domains: converts Array[Struct] to string of comma-separated domain names
        hgvs: set to hgvsp is it exists, or else hgvsc. formats hgvsp for synonymous variants.
        major_consequence: set to most severe consequence for that transcript (
            VEP sometimes provides multiple consequences for a single transcript)
        major_consequence_rank: major_consequence rank based on VEP SO ontology (most severe = 1)
            (see http://www.ensembl.org/info/genome/variation/predicted_data.html)
        category: set to one of: "lof", "missense", "synonymous", "other" based on the value of major_consequence.
    Args:
        vep_root (StructExpression): root path of the VEP struct in the MT
        include_coding_annotations (bool): if True, fields relevant to protein-coding variants will be included
    """

    selected_annotations = [
        "biotype",
        "canonical",
        "cdna_start",
        "cdna_end",
        "codons",
        "gene_id",
        "gene_symbol",
        "hgvsc",
        "hgvsp",
        "transcript_id",
    ]

    if include_coding_annotations:
        selected_annotations.extend(
            [
                "amino_acids",
                "lof",
                "lof_filter",
                "lof_flags",
                "lof_info",
                "polyphen_prediction",
                "protein_id",
                "protein_start",
                "sift_prediction",
            ]
        )

    omit_consequence_terms = (
        hl.set(omit_consequences) if omit_consequences else hl.empty_set(hl.tstr)
    )

    result = hl.sorted(
        vep_root.transcript_consequences.map(
            lambda c: c.select(
                *selected_annotations,
                consequence_terms=c.consequence_terms.filter(
                    lambda t: ~omit_consequence_terms.contains(t)
                ),
                domains=c.domains.map(lambda domain: domain.db + ":" + domain.name),
                major_consequence=hl.cond(
                    c.consequence_terms.size() > 0,
                    hl.sorted(
                        c.consequence_terms,
                        key=lambda t: CONSEQUENCE_TERM_RANK_LOOKUP.get(t),
                    )[0],
                    hl.null(hl.tstr),
                ),
            )
        )
        .filter(lambda c: c.consequence_terms.size() > 0)
        .map(
            lambda c: c.annotate(
                category=(
                    hl.case()
                    .when(
                        CONSEQUENCE_TERM_RANK_LOOKUP.get(c.major_consequence)
                        <= CONSEQUENCE_TERM_RANK_LOOKUP.get("frameshift_variant"),
                        "lof",
                    )
                    .when(
                        CONSEQUENCE_TERM_RANK_LOOKUP.get(c.major_consequence)
                        <= CONSEQUENCE_TERM_RANK_LOOKUP.get("missense_variant"),
                        "missense",
                    )
                    .when(
                        CONSEQUENCE_TERM_RANK_LOOKUP.get(c.major_consequence)
                        <= CONSEQUENCE_TERM_RANK_LOOKUP.get("synonymous_variant"),
                        "synonymous",
                    )
                    .default("other")
                ),
                hgvs=get_expr_for_formatted_hgvs(c),
                major_consequence_rank=CONSEQUENCE_TERM_RANK_LOOKUP.get(
                    c.major_consequence
                ),
            )
        ),
        lambda c: (
            hl.bind(
                lambda is_coding, is_most_severe, is_canonical: (
                    hl.cond(
                        is_coding,
                        hl.cond(
                            is_most_severe,
                            hl.cond(is_canonical, 1, 2),
                            hl.cond(is_canonical, 3, 4),
                        ),
                        hl.cond(
                            is_most_severe,
                            hl.cond(is_canonical, 5, 6),
                            hl.cond(is_canonical, 7, 8),
                        ),
                    )
                ),
                hl.or_else(c.biotype, "") == "protein_coding",
                hl.set(c.consequence_terms).contains(vep_root.most_severe_consequence),
                hl.or_else(c.canonical, 0) == 1,
            )
        ),
    )

    if not include_coding_annotations:
        # for non-coding variants, drop fields here that are hard to exclude in the above code
        result = result.map(lambda c: c.drop("domains", "hgvsp"))

    return hl.zip_with_index(result).map(
        lambda csq_with_index: csq_with_index[1].annotate(
            transcript_rank=csq_with_index[0]
        )
    )


def get_expr_for_vep_protein_domains_set_from_sorted(
    vep_sorted_transcript_consequences_root,
):
    return hl.set(vep_sorted_transcript_consequences_root.flatmap(lambda c: c.domains))


def get_expr_for_vep_gene_id_to_consequence_map(
    vep_sorted_transcript_consequences_root, gene_ids
):
    # Manually build string because hl.json encodes a dictionary as [{ key: ..., value: ... }, ...]
    return (
        "{"
        + hl.delimit(
            gene_ids.map(
                lambda gene_id: hl.bind(
                    lambda worst_consequence_in_gene: '"'
                    + gene_id
                    + '":"'
                    + worst_consequence_in_gene.major_consequence
                    + '"',
                    vep_sorted_transcript_consequences_root.find(
                        lambda c: c.gene_id == gene_id
                    ),
                )
            )
        )
        + "}"
    )


def get_expr_for_vep_transcript_id_to_consequence_map(vep_transcript_consequences_root):
    # Manually build string because hl.json encodes a dictionary as [{ key: ..., value: ... }, ...]
    return (
        "{"
        + hl.delimit(
            vep_transcript_consequences_root.map(
                lambda c: '"' + c.transcript_id + '": "' + c.major_consequence + '"'
            )
        )
        + "}"
    )


def get_expr_for_vep_transcript_ids_set(vep_transcript_consequences_root):
    return hl.set(vep_transcript_consequences_root.map(lambda c: c.transcript_id))


def get_expr_for_worst_transcript_consequence_annotations_struct(
    vep_sorted_transcript_consequences_root, include_coding_annotations=True
):
    """Retrieves the top-ranked transcript annotation based on the ranking computed by
    get_expr_for_vep_sorted_transcript_consequences_array(..)
    Args:
        vep_sorted_transcript_consequences_root (ArrayExpression):
        include_coding_annotations (bool):
    """

    transcript_consequences = {
        "biotype": hl.tstr,
        "canonical": hl.tint,
        "category": hl.tstr,
        "cdna_start": hl.tint,
        "cdna_end": hl.tint,
        "codons": hl.tstr,
        "gene_id": hl.tstr,
        "gene_symbol": hl.tstr,
        "hgvs": hl.tstr,
        "hgvsc": hl.tstr,
        "major_consequence": hl.tstr,
        "major_consequence_rank": hl.tint,
        "transcript_id": hl.tstr,
    }

    if include_coding_annotations:
        transcript_consequences.update(
            {
                "amino_acids": hl.tstr,
                "domains": hl.tstr,
                "hgvsp": hl.tstr,
                "lof": hl.tstr,
                "lof_flags": hl.tstr,
                "lof_filter": hl.tstr,
                "lof_info": hl.tstr,
                "polyphen_prediction": hl.tstr,
                "protein_id": hl.tstr,
                "sift_prediction": hl.tstr,
            }
        )

    return hl.cond(
        vep_sorted_transcript_consequences_root.size() == 0,
        hl.struct(
            **{
                field: hl.null(field_type)
                for field, field_type in transcript_consequences.items()
            }
        ),
        hl.bind(
            lambda worst_transcript_consequence: (
                worst_transcript_consequence.annotate(
                    domains=hl.delimit(hl.set(worst_transcript_consequence.domains))
                ).select(*transcript_consequences.keys())
            ),
            vep_sorted_transcript_consequences_root[0],
        ),
    )


def get_expr_for_alt_allele(table: hl.Table) -> hl.str:
    return table.alleles[1]


def get_expr_for_contig(locus: hl.expr.LocusExpression) -> hl.expr.StringExpression:
    """Normalized contig name"""
    return locus.contig.replace("^chr", "")


def get_expr_for_contig_number(
    locus: hl.expr.LocusExpression,
) -> hl.expr.Int32Expression:
    """Convert contig name to contig number"""
    return hl.bind(
        lambda contig: (
            hl.case()
            .when(contig == "X", 23)
            .when(contig == "Y", 24)
            .when(contig[0] == "M", 25)
            .default(hl.int(contig))
        ),
        get_expr_for_contig(locus),
    )


def get_expr_for_variant_ids(
    locus: hl.expr.LocusExpression,
    alleles: hl.expr.ArrayExpression,
    max_length: int = None,
) -> hl.expr.ArrayExpression:
    """Return a list of variant ids - one for each alt allele in the variant"""

    def compute_variant_id(alt):
        variant_id = (
            locus.contig + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alt
        )
        if max_length is not None:
            variant_id = variant_id[:max_length]
        return variant_id

    return alleles[1:].map(compute_variant_id)


def get_expr_for_variant_type(table: hl.Table) -> hl.str:
    return hl.bind(
        lambda ref_len, alt_len: (
            hl.case()
            .when(ref_len > alt_len, "D")
            .when(ref_len < alt_len, "I")
            .when(ref_len > 1, "M")
            .default("S")
        ),
        hl.len(get_expr_for_ref_allele(table)),
        hl.len(get_expr_for_alt_allele(table)),
    )


def get_expr_for_ref_allele(table):
    return table.alleles[0]


def get_expr_for_start_pos(table):
    return table.locus.position


def get_expr_for_end_pos(table):
    return table.locus.position + hl.len(get_expr_for_ref_allele(table)) - 1


def get_expr_for_variant_id(table, max_length=None):
    """Expression for computing <chrom>-<pos>-<ref>-<alt>. Assumes alleles were split.
    Args:
        max_length: (optional) length at which to truncate the <chrom>-<pos>-<ref>-<alt> string
    Return:
        string: "<chrom>-<pos>-<ref>-<alt>"
    """
    contig = get_expr_for_contig(table.locus)
    variant_id = (
        contig
        + "-"
        + hl.str(table.locus.position)
        + "-"
        + table.alleles[0]
        + "-"
        + table.alleles[1]
    )
    if max_length is not None:
        return variant_id[0:max_length]
    return variant_id


def get_expr_for_xpos(locus: hl.expr.LocusExpression) -> hl.expr.Int64Expression:
    """Genomic position represented as a single number = contig_number * 10**9 + position.
    This represents chrom:pos more compactly and allows for easier sorting.
    """
    contig_number = get_expr_for_contig_number(locus)
    return hl.int64(contig_number) * 1_000_000_000 + locus.position


@click.command()
@click.option(
    '--vcf-path',
    'vcf_path',
    required=True,
)
@click.option(
    '--siteonly-vqsr-vcf-path',
    'siteonly_vqsr_vcf_path',
    required=True,
)
@click.option(
    '--vep-ht-path',
    'vep_ht_path',
    required=True,
)
@click.option(
    '--out-mt-path',
    'out_mt_path',
    required=True,
)
@click.option(
    '--checkpoint-prefix',
    'checkpoint_prefix',
    required=True,
)
def main(
    vcf_path: str,
    siteonly_vqsr_vcf_path: str,
    vep_ht_path: str,
    out_mt_path: str,
    checkpoint_prefix: str,
):
    hl.init(default_reference='GRCh38')

    annotate_cohort(
        vcf_path=vcf_path,
        site_only_vqsr_vcf_path=siteonly_vqsr_vcf_path,
        vep_ht_path=vep_ht_path,
        out_mt_path=out_mt_path,
        overwrite=not get_config()['workflow'].get('check_intermediates'),
        sequencing_type=get_config()['workflow']['sequencing_type'],
        checkpoint_prefix=checkpoint_prefix,
    )


def load_vqsr(site_only_vqsr_vcf_path: str, genome_build: str = 'GRCh38'):
    """
    Loads the site-only VCF with applied AS-VQSR filters into a site-only hail table.
    Populates "ht.filters".
    """
    logging.info(
        f'AS-VQSR: importing annotations from a site-only VCF '
        f'{site_only_vqsr_vcf_path}'
    )
    ht = hl.import_vcf(
        str(site_only_vqsr_vcf_path),
        force_bgz=True,
        reference_genome=genome_build,
    ).rows()

    # Annotating FILTERS from AS_FilterStatus before dropping that latter
    ht = ht.annotate(
        filters=ht.filters.union(hl.set(ht.info.AS_FilterStatus)).filter(
            lambda val: val != 'PASS'
        ),
    )

    # Dropping all INFO/AS* annotations as they are causing problems parsing by Hail.
    as_fields = [f for f in ht.info if f.startswith('AS_')]
    as_fields_to_keep: list[str] = []
    as_fields_to_drop = [f for f in as_fields if f not in as_fields_to_keep]
    ht = ht.annotate(info=ht.info.drop(*as_fields_to_drop))

    logging.info(f'AS-VQSR: splitting multiallelics...')
    unsplit_count = ht.count()
    ht = hl.split_multi_hts(ht)

    split_count = ht.count()
    logging.info(
        f'AS-VQSR: Found {unsplit_count} unsplit and {split_count} split variants'
    )
    return ht


def annotate_cohort(
    vcf_path,
    site_only_vqsr_vcf_path,
    vep_ht_path,
    out_mt_path,
    overwrite=False,
    sequencing_type='',
    checkpoint_prefix=None,
):
    """
    Convert VCF to mt, annotate for seqr loader, add VEP annotations.
    """

    def _read(path):
        if path.strip('/').endswith('.ht'):
            t = hl.read_table(str(path))
        else:
            assert path.strip('/').endswith('.mt')
            t = hl.read_matrix_table(str(path))
        logging.info(f'Read checkpoint {path}')
        return t

    def _checkpoint(t, file_name):
        if checkpoint_prefix:
            path = os.path.join(checkpoint_prefix, file_name)
            if not overwrite and hl.hadoop_exists(os.path.join(path, '_SUCCESS')):
                t = _read(str(path))
            else:
                t.write(str(path), overwrite=True)
                logging.info(f'Wrote checkpoint {path}')
                t = _read(str(path))
        return t

    mt = hl.import_vcf(
        str(vcf_path),
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
    )
    logging.info(
        f'Importing VCF {vcf_path}, '
        f'adding VQSR annotations from {site_only_vqsr_vcf_path}, '
        f'adding VEP annotations from {vep_ht_path}'
    )

    logging.info(f'Loading VEP Table from {vep_ht_path}')
    # Annotate VEP. Do ti before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(str(vep_ht_path))
    logging.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # Splitting multi-allelics. We do not handle AS info fields here - we handle
    # them when loading VQSR instead, and populate entrie "info" from VQSR.
    mt = hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )
    mt = _checkpoint(mt, 'mt-vep-split.mt')

    vqsr_ht = load_vqsr(site_only_vqsr_vcf_path)
    overwrite = True
    vqsr_ht = _checkpoint(vqsr_ht, 'vqsr.ht')

    logging.info('Adding VQSR annotations into the Matrix Table')
    mt = mt.annotate_globals(**vqsr_ht.index_globals())
    mt = mt.annotate_rows(
        # vqsr_ht has info annotation split by allele, plus the new AS-VQSR annotations
        info=vqsr_ht[mt.row_key].info,
        filters=mt.filters.union(vqsr_ht[mt.row_key].filters).filter(
            lambda val: val != 'PASS'
        ),
    )
    mt = _checkpoint(mt, 'mt-vep-split-vqsr.mt')
    overwrite = False

    ref_ht = hl.read_table(str(reference_path('seqr/combined_reference')))
    clinvar_ht = hl.read_table(str(reference_path('seqr/clinvar')))

    logging.info('Annotating with seqr-loader fields: round 1')
    mt = mt.annotate_rows(
        AC=mt.info.AC[mt.a_index - 1],
        AF=mt.info.AF[mt.a_index - 1],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=get_expr_for_variant_ids(mt.locus_old, mt.alleles_old),
        sortedTranscriptConsequences=get_expr_for_vep_sorted_transcript_consequences_array(
            mt.vep
        ),
        variantId=get_expr_for_variant_id(mt),
        contig=get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
    )
    mt = _checkpoint(mt, 'mt-vep-split-vqsr-round1.mt')

    logging.info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)'
    )
    mt = mt.annotate_rows(
        domains=get_expr_for_vep_protein_domains_set_from_sorted(
            mt.sortedTranscriptConsequences
        ),
        transcriptConsequenceTerms=get_expr_for_vep_consequence_terms_set(
            mt.sortedTranscriptConsequences
        ),
        transcriptIds=get_expr_for_vep_transcript_ids_set(
            mt.sortedTranscriptConsequences
        ),
        mainTranscript=get_expr_for_worst_transcript_consequence_annotations_struct(
            mt.sortedTranscriptConsequences
        ),
        geneIds=get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences),
        codingGeneIds=get_expr_for_vep_gene_ids_set(
            mt.sortedTranscriptConsequences, only_coding_genes=True
        ),
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
            }
        ),
    )
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
    )
    if sequencing_type:
        # Map to Seqr-style string
        # https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        mt = mt.annotate_globals(
            sampleType={
                'genome': 'WGS',
                'exome': 'WES',
                'single_cell': 'RNA',
            }.get(sequencing_type, ''),
        )

    logging.info('Done:')
    mt.describe()
    mt.write(str(out_mt_path), overwrite=overwrite)
    logging.info(f'Written final matrix table into {out_mt_path}')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
