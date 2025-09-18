import itertools
import logging

import hail as hl

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import genome_build
from cpg_workflows.utils import can_reuse
from gnomad.utils.filtering import add_filters_expr

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

EMPTY_TABLE_CONFIG = """
array<
    struct{
        locus: locus<GRCh38>,
        alleles: array<str>,
        variant_id: str,
        rsids: str,
        type: struct {
            colocated_variants: struct {
                all: array<str>
            },
            subsets: set<str>,
            flags: set<str>,
            freq: struct {
                all: struct {
                    ac: int32,
                    an: int32,
                    hemizygote_count: int32,
                    homozygote_count: int32,
                    ancestry_groups: array<struct {
                        id: str,
                        ac: int32,
                        an: int32,
                        hemizygote_count: int32,
                        homozygote_count: int32
                    }>
                }
            },
            fafmax: struct {
                faf95_max: float64,
                faf95_max_gen_anc: str,
                faf99_max: float64,
                faf99_max_gen_anc: str
            },
            vep: struct {
                minimised: int32,
                assembly_name: str,
                allele_string: str,
                ancestral: str,
                colocated_variants: array<struct {
                    allele_string: str,
                    clin_sig: array<str>,
                    clin_sig_allele: str,
                    end: int32,
                    id: str,
                    minimised: int32,
                    minor_allele: str,
                    minor_allele_freq: float64,
                    phenotype_or_disease: int32,
                    pubmed: array<int32>,
                    seq_region_name: str,
                    somatic: int32,
                    start: int32,
                    strand: int32
                }>,
                context: str,
                end: int32,
                id: str,
                input: str,
                intergenic_consequences: array<struct {
                    allele_num: int32,
                    consequence_terms: array<str>,
                    impact: str,
                    minimised: int32,
                    variant_allele: str
                }>,
                most_severe_consequence: str,
                motif_feature_consequences: array<struct {
                    allele_num: int32,
                    consequence_terms: array<str>,
                    high_inf_pos: str,
                    impact: str,
                    minimised: int32,
                    motif_feature_id: str,
                    motif_name: str,
                    motif_pos: int32,
                    motif_score_change: float64,
                    strand: int32,
                    transcription_factors: array<str>,
                    variant_allele: str
                }>,
                regulatory_feature_consequences: array<struct {
                    allele_num: int32,
                    biotype: str,
                    consequence_terms: array<str>,
                    impact: str,
                    minimised: int32,
                    regulatory_feature_id: str,
                    variant_allele: str
                }>,
                seq_region_name: str,
                start: int32,
                strand: int32,
                transcript_consequences: array<struct {
                    allele_num: int32,
                    amino_acids: str,
                    appris: str,
                    biotype: str,
                    canonical: int32,
                    mane_select: str,
                    mane_plus_clinical: str,
                    ccds: str,
                    cdna_start: int32,
                    cdna_end: int32,
                    cds_end: int32,
                    cds_start: int32,
                    codons: str,
                    consequence_terms: array<str>,
                    distance: int32,
                    domains: array<struct {
                        db: str,
                        name: str
                    }>,
                    exon: str,
                    gene_id: str,
                    gene_pheno: int32,
                    gene_symbol: str,
                    gene_symbol_source: str,
                    hgnc_id: str,
                    hgvsc: str,
                    hgvsp: str,
                    hgvs_offset: int32,
                    impact: str,
                    intron: str,
                    lof: str,
                    lof_flags: str,
                    lof_filter: str,
                    lof_info: str,
                    existing_inframe_oorfs: int32,
                    existing_outofframe_oorfs: int32,
                    existing_uorfs: int32,
                    5utr_consequence: str,
                    5utr_annotation: dict<str, struct {
                        type: str,
                        KozakContext: str,
                        KozakStrength: str,
                        DistanceToCDS: str,
                        CapDistanceToStart: str,
                        DistanceToStop: str,
                        Evidence: str,
                        AltStop: str,
                        AltStopDistanceToCDS: str,
                        FrameWithCDS: str,
                        StartDistanceToCDS: str,
                        newSTOPDistanceToCDS: str,
                        alt_type: str,
                        alt_type_length: str,
                        ref_StartDistanceToCDS: str,
                        ref_type: str,
                        ref_type_length: str
                    }>,
                    minimised: int32,
                    mirna: array<str>,
                    polyphen_prediction: str,
                    polyphen_score: float64,
                    protein_end: int32,
                    protein_start: int32,
                    protein_id: str,
                    sift_prediction: str,
                    sift_score: float64,
                    strand: int32,
                    swissprot: array<str>,
                    transcript_id: str,
                    trembl: array<str>,
                    tsl: int32,
                    uniparc: array<str>,
                    uniprot_isoform: array<str>,
                    variant_allele: str,
                    am_class: str,
                    am_pathogenicity: float64,
                    source: str,
                    flags: array<str>
                }>,
                variant_class: str
            },
            age_distribution: struct {
                het: struct {
                    bin_edges: array<float64>,
                    bin_freq: array<int64>,
                    n_smaller: int64,
                    n_larger: int64
                },
                hom: struct {
                    bin_edges: array<float64>,
                    bin_freq: array<int64>,
                    n_smaller: int64,
                    n_larger: int64
                }
            },
            filters: set<str>,
            quality_metrics: struct {
                allele_balance: struct {
                    alt_adj: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    alt_raw: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    }
                },
                genotype_depth: struct {
                    all_adj: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    all_raw: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    alt_adj: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    alt_raw: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    }
                },
                genotype_quality: struct {
                    all_adj: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    all_raw: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    alt_adj: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    },
                    alt_raw: struct {
                        bin_edges: array<float64>,
                        bin_freq: array<int64>,
                        n_smaller: int64,
                        n_larger: int64
                    }
                },
                site_quality_metrics: array<struct {
                    metric: str,
                    value: float64
                }>
            }
        }
    }
>
"""

OMIT_CONSEQUENCE_TERMS = hl.set(["upstream_gene_variant", "downstream_gene_variant"])

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
    },
)

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


# hail DictExpression that maps each CONSEQUENCE_TERM to its rank in the list
CONSEQUENCE_TERM_RANK_LOOKUP = hl.dict({term: rank for rank, term in enumerate(CONSEQUENCE_TERMS)})


def _nullify_nan(value):
    return hl.if_else(hl.is_nan(value), hl.missing(value.dtype), value)


def _normalized_contig(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.rbind(hl.str(contig).replace("^chr", ""), lambda c: hl.if_else(c == "MT", "M", c))


def _variant_id(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int | None = None):
    """
    Expression for computing <chrom>-<pos>-<ref>-<alt>. Assumes alleles were split.

    Args:
        max_length: (optional) length at which to truncate the <chrom>-<pos>-<ref>-<alt> string

    Return:
        string: "<chrom>-<pos>-<ref>-<alt>"
    """
    contig = _normalized_contig(locus.contig)
    var_id = contig + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alleles[1]

    if max_length is not None:
        return var_id[0:max_length]

    return var_id


def _freq_index_key(subset=None, pop=None, sex=None, raw=False):
    parts = [s for s in [subset, pop, sex] if s is not None]
    parts.append("raw" if raw else "adj")
    return "_".join(parts)


def _freq(ds, *args, **kwargs):
    return ds.freq[ds.freq_index_dict[_freq_index_key(*args, **kwargs)]]


def _subset_filter(subset):
    return lambda variant: variant.ac_adj[subset] > 0


def process_score_cutoffs(
    snv_bin_cutoff: int | None = None,
    indel_bin_cutoff: int | None = None,
    snv_score_cutoff: float | None = None,
    indel_score_cutoff: float | None = None,
    aggregated_bin_ht_path: hl.Table = None,
    snv_bin_id: str = 'bin',
    indel_bin_id: str = 'bin',
) -> dict[str, hl.expr.StructExpression]:
    """
    NOTE: This function was lifted from the gnomad_qc repo `gnomad_qc.v4.variant_qc.final_filter.process_score_cutoffs`.
        - Changes made:
            - Removed `ht` parameter, as it is not needed.
                - Global `min_score` and `max_score` are now calculated from the `aggregated_bin_ht`
    Determine SNP and indel score cutoffs if given bin instead of score.

    .. note::

        - `snv_bin_cutoff` and `snv_score_cutoff` are mutually exclusive, and one must
          be supplied.
        - `indel_bin_cutoff` and `indel_score_cutoff` are mutually exclusive, and one
          must be supplied.
        - If a `snv_bin_cutoff` or `indel_bin_cutoff` cutoff is supplied then an
          `aggregated_bin_ht` and `bin_id` must also be supplied to determine the SNP
          and indel scores to use as cutoffs from an aggregated bin Table like one
          created by `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param ht: Filtering Table to prepare as the final filter Table.
    :param snv_bin_cutoff: Bin cutoff to use for SNP variant QC filter. Can't be used
        with `snv_score_cutoff`.
    :param indel_bin_cutoff: Bin cutoff to use for indel variant QC filter. Can't be
        used with `indel_score_cutoff`.
    :param snv_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for SNP variant QC filter. Can't be used with `snv_bin_cutoff`.
    :param indel_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for indel variant QC filter. Can't be used with `indel_bin_cutoff`.
    :param aggregated_bin_ht: Table with aggregate counts of variants based on bins
    :param snv_bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to
        determine the SNP score cutoff.
    :param indel_bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to
        determine the indel score cutoff.
    :return: Finalized random forest Table annotated with variant filters.
    """
    aggregated_bin_ht = hl.read_table(aggregated_bin_ht_path)

    if snv_bin_cutoff is not None and snv_score_cutoff is not None:
        raise ValueError(
            "snv_bin_cutoff and snv_score_cutoff are mutually exclusive, please only"
            " supply one SNP filtering cutoff.",
        )

    if indel_bin_cutoff is not None and indel_score_cutoff is not None:
        raise ValueError(
            "indel_bin_cutoff and indel_score_cutoff are mutually exclusive, please"
            " only supply one indel filtering cutoff.",
        )

    if snv_bin_cutoff is None and snv_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters snv_bin_cutoff and snv_score_cutoff" " must be supplied.",
        )

    if indel_bin_cutoff is None and indel_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters indel_bin_cutoff and" " indel_score_cutoff must be supplied.",
        )

    if (
        (snv_bin_cutoff is not None and snv_bin_id is None) or (indel_bin_cutoff is not None and indel_bin_id is None)
    ) and (aggregated_bin_ht is None):
        raise ValueError(
            "If using snv_bin_cutoff or indel_bin_cutoff, both aggregated_bin_ht and"
            " snv_bin_id/indel_bin_id must be supplied",
        )

    cutoffs = {
        "snv": {"score": snv_score_cutoff, "bin": snv_bin_cutoff, "bin_id": snv_bin_id},
        "indel": {
            "score": indel_score_cutoff,
            "bin": indel_bin_cutoff,
            "bin_id": indel_bin_id,
        },
    }

    # Calculate min and max scores within each bin.
    min_score = aggregated_bin_ht.aggregate(hl.agg.min(aggregated_bin_ht.min_score))
    max_score = aggregated_bin_ht.aggregate(hl.agg.max(aggregated_bin_ht.max_score))
    aggregated_bin_ht = aggregated_bin_ht.annotate(indel=~aggregated_bin_ht.snv)

    cutoff_globals = {}
    for variant_type, cutoff in cutoffs.items():
        if cutoff["bin"] is not None:
            score_cutoff = aggregated_bin_ht.aggregate(
                hl.agg.filter(
                    aggregated_bin_ht[variant_type]
                    & (aggregated_bin_ht.bin_id == cutoff["bin_id"])
                    & (aggregated_bin_ht.bin == cutoff["bin"]),
                    hl.agg.min(aggregated_bin_ht.min_score),
                ),
            )
            cutoff_globals[variant_type] = hl.struct(
                bin=cutoff["bin"],
                min_score=score_cutoff,
                bin_id=cutoff["bin_id"],
            )
        else:
            cutoff_globals[variant_type] = hl.struct(min_score=cutoff["score"])

        score_cutoff = hl.eval(cutoff_globals[variant_type].min_score)
        if score_cutoff < min_score or score_cutoff > max_score:
            raise ValueError(
                f"{variant_type}_score_cutoff is not within the range of score (" f"{min_score, max_score}).",
            )

    logger.info(
        f"Using a SNP score cutoff of {hl.eval(cutoff_globals['snv'].min_score)} and an"
        f" indel score cutoff of {hl.eval(cutoff_globals['indel'].min_score)}.",
    )

    return cutoff_globals


def consequence_term_rank(consequence_term):
    return CONSEQUENCE_TERM_RANK_LOOKUP.get(consequence_term)


def hgvsp_from_consequence_amino_acids(csq):
    return hl.if_else(
        csq.hgvsp.contains("=") | csq.hgvsp.contains("%3D"),
        hl.bind(
            lambda protein_letters: "p." + protein_letters + hl.str(csq.protein_start) + protein_letters,
            hl.delimit(
                csq.amino_acids.split("").filter(lambda letter: letter != "").map(PROTEIN_LETTERS_1TO3.get),
                "",
            ),
        ),
        csq.hgvsp.split(":")[-1],
    )


def extract_transcripts(ds: hl.Table) -> hl.Table:
    ds = ds.key_by()
    ds = ds.select(gene=ds.row_value)
    ds = ds.annotate(transcripts=ds.gene.transcripts)
    ds = ds.explode(ds.transcripts)
    ds = ds.annotate(**ds.transcripts).drop("transcripts")
    ds = ds.key_by("transcript_id")
    ds = ds.repartition(2000, shuffle=True)

    return ds


def annotate_transcript_consequences(ds: hl.Table, transcripts: hl.Table, mane_transcripts_path=None):

    most_severe_consequence = ds.vep.most_severe_consequence

    transcript_consequences = ds.vep.transcript_consequences

    # Drop irrelevant consequences
    transcript_consequences = transcript_consequences.map(
        lambda c: c.annotate(
            consequence_terms=c.consequence_terms.filter(
                lambda t: ~OMIT_CONSEQUENCE_TERMS.contains(t),  # pylint: disable=invalid-unary-operand-type
            ),
        ),
    ).filter(lambda c: c.consequence_terms.size() > 0)

    # Add/transmute derived fields
    transcript_consequences = transcript_consequences.map(
        lambda c: c.annotate(major_consequence=hl.sorted(c.consequence_terms, key=consequence_term_rank)[0]),
    ).map(
        lambda c: c.annotate(
            domains=hl.set(c.domains.map(lambda domain: domain.db + ":" + domain.name).filter(hl.is_defined)),
            hgvsc=c.hgvsc.split(":")[-1],
            hgvsp=hgvsp_from_consequence_amino_acids(c),
            is_canonical=hl.bool(c.canonical),
        ),
    )

    consequences = [
        "biotype",
        "consequence_terms",
        "domains",
        "gene_id",
        "gene_symbol",
        "hgvsc",
        "hgvsp",
        "is_canonical",
        "lof_filter",
        "lof_flags",
        "lof",
        "major_consequence",
        "transcript_id",
    ]

    # gnomAD v2 and v3 have these consequences. gnomAD v4 does not.
    v2_and_v3_only_consequences = ["polyphen_prediction", "sift_prediction"]

    available_consequences = list(transcript_consequences[0])

    for csq in v2_and_v3_only_consequences:
        if csq in available_consequences:
            consequences.append(csq)

    transcript_consequences = transcript_consequences.map(lambda c: c.select(*consequences))

    # TODO: This can potentially be improved by removing Table.collect
    # See https://hail.zulipchat.com/#narrow/stream/123010-Hail-0.2E2.20support/topic/Optimize.20annotation.20with.20small.20dataset
    # and https://github.com/Nealelab/ukb_common/blob/ad94d20f8c9f3b711e40a473425925775f0b1f30/utils/generic.py#L18
    transcript_info = hl.dict(
        [
            (row.transcript_id, row.transcript_info)
            for row in transcripts.select(
                transcript_info=hl.struct(
                    transcript_version=transcripts.transcript_version,
                    gene_version=transcripts.gene.gene_version,
                ),
            ).collect()
        ],
    )

    transcript_consequences = transcript_consequences.map(
        lambda csq: csq.annotate(**transcript_info.get(csq.transcript_id)),
    )

    if mane_transcripts_path:
        mane_transcripts = hl.read_table(mane_transcripts_path)
        mane_select_transcripts_version = hl.eval(mane_transcripts.globals.version)

        mane_transcripts = hl.dict([(row.gene_id, row.drop("gene_id")) for row in mane_transcripts.collect()])

        transcript_consequences = transcript_consequences.map(
            lambda csq: csq.annotate(
                **hl.rbind(
                    mane_transcripts.get(csq.gene_id),
                    lambda mane_transcript: (
                        hl.case()
                        .when(
                            (mane_transcript.ensembl_id == csq.transcript_id)
                            & (mane_transcript.ensembl_version == csq.transcript_version),
                            hl.struct(
                                is_mane_select=True,
                                is_mane_select_version=True,
                                refseq_id=mane_transcript.refseq_id,
                                refseq_version=mane_transcript.refseq_version,
                            ),
                        )
                        .when(
                            mane_transcript.ensembl_id == csq.transcript_id,
                            hl.struct(
                                is_mane_select=True,
                                is_mane_select_version=False,
                                refseq_id=hl.null(hl.tstr),
                                refseq_version=hl.null(hl.tstr),
                            ),
                        )
                        .default(
                            hl.struct(
                                is_mane_select=False,
                                is_mane_select_version=False,
                                refseq_id=hl.null(hl.tstr),
                                refseq_version=hl.null(hl.tstr),
                            ),
                        )
                    ),
                ),
            ),
        )

        transcript_consequences = hl.sorted(
            transcript_consequences,
            lambda c: (
                hl.if_else(c.biotype == "protein_coding", 0, 1, missing_false=True),
                hl.if_else(c.major_consequence == most_severe_consequence, 0, 1, missing_false=True),
                hl.if_else(c.is_mane_select, 0, 1, missing_false=True),
                hl.if_else(c.is_canonical, 0, 1, missing_false=True),
            ),
        )

        ds = ds.annotate(transcript_consequences=transcript_consequences).drop("vep")
        ds = ds.annotate_globals(mane_select_version=mane_select_transcripts_version)

    else:
        transcript_consequences = hl.sorted(
            transcript_consequences,
            lambda c: (
                hl.if_else(c.biotype == "protein_coding", 0, 1, missing_false=True),
                hl.if_else(c.major_consequence == most_severe_consequence, 0, 1, missing_false=True),
                hl.if_else(c.is_canonical, 0, 1, missing_false=True),
            ),
        )

        ds = ds.annotate(transcript_consequences=transcript_consequences).drop("vep")

    return ds


def adjust_interval_padding(ht: hl.Table, padding: int) -> hl.Table:
    """
    Adjust interval padding in HT.

    .. warning::

        This function can lead to overlapping intervals, so it is not recommended for
        most applications. For example, it can be used to filter a variant list to all
        variants within the returned interval list, but would not work for getting an
        aggregate statistic for each interval if the desired output is independent
        statistics.

    :param ht: HT to adjust.
    :param padding: Padding to use.
    :return: HT with adjusted interval padding.
    """
    return ht.key_by(
        interval=hl.locus_interval(
            ht.interval.start.contig,
            ht.interval.start.position - padding,
            ht.interval.end.position + padding,
            # Include the end of the intervals to capture all variants.
            reference_genome=ht.interval.start.dtype.reference_genome,
            includes_end=True,
        ),
    )


def prepare_gnomad_v4_variants_helper(
    ds_path: str | None,
    exome_or_genome: str,
    score_cutoffs: dict[str, hl.expr.StructExpression],
) -> hl.Table:

    if ds_path is None:
        ds = hl.Table.parallelize(hl.literal([], "".join(EMPTY_TABLE_CONFIG.split()))).key_by('locus', 'alleles')

        age_hist = hl.struct(
            bin_edges=hl.literal([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0]),
            bin_freq=hl.literal([0] * 10),
            n_smaller=hl.int32(0),
            n_larger=hl.int32(0),
        )

        ds_globals = ds.globals.annotate(**{f"{exome_or_genome}_age_distribution": age_hist})
        ds = ds.annotate_globals(**ds_globals)

        ds = ds.rename({'type': exome_or_genome})
        return ds

    logger.info('Preparing the base frequencies table...')
    ds = hl.read_table(ds_path)

    subsets = set(m.get("subset", None) for m in ds.freq_meta.collect()[0])

    ds = ds.annotate(variant_id=_variant_id(ds.locus, ds.alleles))
    ds = ds.rename({"rsid": "rsids"})

    # Identify co-located variants.
    logger.info('Identifying co-located variants...')
    variants_by_locus = ds.select(
        ds.variant_id,
        ac_adj=hl.struct(**{subset or "all": _freq(ds, subset=subset).AC for subset in subsets}),
    )

    variants_by_locus = variants_by_locus.group_by("locus").aggregate(
        variants=hl.agg.collect(variants_by_locus.row_value),
    )

    variants_by_locus = variants_by_locus.annotate(
        variant_ids=hl.struct(
            **{
                subset
                or "all": variants_by_locus.variants.filter(_subset_filter(subset or "all")).map(
                    lambda variant: variant.variant_id,
                )
                for subset in subsets
            },
        ),
    )

    ds = ds.annotate(colocated_variants=variants_by_locus[ds.locus].variant_ids)
    ds = ds.annotate(
        colocated_variants=hl.struct(
            **{
                subset: ds.colocated_variants[subset].filter(lambda variant_id: variant_id != ds.variant_id)
                for subset in ds.colocated_variants._fields
            },
        ),
    )

    # Expand the frequencies struct.
    logger.info('Expanding the frequencies struct...')

    # Annotate when a locus is not on the autosomes or in the PAR region, and
    # can therefore have hemizygote calls for XY individuals.
    ds = ds.annotate(in_autosome_or_par=ds.locus.in_autosome_or_par())

    # Identify the ancestry groups in each subset.
    subset_ancestry_groups = {}
    for subset in subsets:
        subset_ancestry_groups[subset] = set(
            m.get("pop", None) for m in ds.freq_meta.collect()[0] if m.get("subset", None) == subset
        )
        subset_ancestry_groups[subset].discard(None)

    # Annotate the dataset with a frequency struct describing the call stats for each
    # subset, optionally split by ancestry group and sex.
    ds = ds.annotate(
        expanded_freq=hl.struct(
            **{
                subset
                or "all": hl.struct(
                    ac=_freq(ds, subset=subset).AC,
                    an=_freq(ds, subset=subset).AN,
                    hemizygote_count=hl.if_else(
                        ds.in_autosome_or_par,
                        0,
                        hl.or_else(_freq(ds, subset=subset, sex="XY").AC, 0),
                    ),
                    homozygote_count=_freq(ds, subset=subset).homozygote_count,
                    ancestry_groups=[
                        hl.struct(
                            id="_".join(filter(bool, [pop, sex])),  # type: ignore
                            ac=hl.or_else(_freq(ds, subset=subset, pop=pop, sex=sex).AC, 0),
                            an=hl.or_else(_freq(ds, subset=subset, pop=pop, sex=sex).AN, 0),
                            hemizygote_count=(
                                0
                                if sex == "XX"
                                else hl.if_else(
                                    ds.in_autosome_or_par,
                                    0,
                                    hl.or_else(_freq(ds, subset=subset, pop=pop, sex="XY").AC, 0),
                                )
                            ),
                            homozygote_count=hl.or_else(
                                _freq(ds, subset=subset, pop=pop, sex=sex).homozygote_count,
                                0,
                            ),
                        )
                        for pop, sex in list(itertools.product(subset_ancestry_groups[subset], [None, "XX", "XY"]))
                        + [(None, "XX"), (None, "XY")]
                    ],
                )
                for subset in subsets
            },
        ),
    )

    # Flag subsets in which the variant is present
    logger.info('Identifying the subsets that contain each variant...')
    filtered_subsets = [subset for subset in subsets if subset is not None]
    if not filtered_subsets:
        ds = ds.annotate(subsets=hl.empty_set(hl.tstr))
    else:
        ds = ds.annotate(
            subsets=hl.set(
                hl.array([(subset, ds.expanded_freq[subset].ac > 0) for subset in filtered_subsets])
                .filter(lambda t: t[1])
                .map(lambda t: t[0]),
            ),
        )

    # Summarise the carrier age distribution for each variant.
    logger.info('Summarising the age distribution for variant carriers...')
    ds = ds.annotate(expanded_age_distribution=hl.struct(het=ds.age_hist_het, hom=ds.age_hist_hom))

    # Summarise the site quality metrics.
    logger.info('Summarising the site quality metrics...')
    ds = ds.annotate(
        quality_metrics=hl.struct(
            allele_balance=hl.struct(
                alt_adj=ds.histograms.qual_hists.ab_hist_alt,
                alt_raw=ds.histograms.raw_qual_hists.ab_hist_alt,
            ),
            genotype_depth=hl.struct(
                all_adj=ds.histograms.qual_hists.dp_hist_all,
                all_raw=ds.histograms.raw_qual_hists.dp_hist_all,
                alt_adj=ds.histograms.qual_hists.dp_hist_alt,
                alt_raw=ds.histograms.raw_qual_hists.dp_hist_alt,
            ),
            genotype_quality=hl.struct(
                all_adj=ds.histograms.qual_hists.gq_hist_all,
                all_raw=ds.histograms.raw_qual_hists.gq_hist_all,
                alt_adj=ds.histograms.qual_hists.gq_hist_alt,
                alt_raw=ds.histograms.raw_qual_hists.gq_hist_alt,
            ),
            site_quality_metrics=[hl.struct(metric="SiteQuality", value=hl.float(_nullify_nan(ds.info.QUALapprox)))]
            + [hl.struct(metric="inbreeding_coeff", value=hl.float(_nullify_nan(ds.inbreeding_coeff[0])))]
            + [
                hl.struct(
                    metric=metric,
                    value=hl.float(
                        _nullify_nan(ds.info[metric]),
                    ),  # default_compute_info() annotates below fields as hl.array(hl.float64) so need to get value in array
                )
                for metric in [
                    "AS_MQ",
                    "AS_FS",
                    "AS_MQRankSum",
                    "AS_pab_max",
                    "AS_QUALapprox",
                    "AS_QD",
                    "AS_ReadPosRankSum",
                    "AS_SOR",
                    "AS_VarDP",
                    "AS_VQSLOD",
                ]
            ],
        ),
    )

    # Restructure the region flags.
    logger.info('Restructuring the region flags...')
    flags = [
        hl.or_missing(ds.region_flags.lcr, "lcr"),
        hl.or_missing(ds.region_flags.segdup, "segdup"),
        hl.or_missing(ds.region_flags.non_par, "non_par"),
        hl.or_missing(
            ((ds.locus.contig == "chrX") & ds.locus.in_x_par()) | ((ds.locus.contig == "chrY") & ds.locus.in_y_par()),
            "par",
        ),
        hl.or_missing(ds.monoallelic, "monoallelic"),
    ]
    ds = ds.annotate(flags=hl.set(flags).filter(hl.is_defined))

    # Add the variant filter flags.
    logger.info('Adding the variant filter flags...')

    inbreeding_coeff_cutoff = config_retrieve(['large_cohort', 'browser', 'inbreeding_coeff_cutoff'])
    filters = {
        "InbreedingCoeff": ds.inbreeding_coeff[0] < inbreeding_coeff_cutoff,
        "AC0": ds.expanded_freq.all.ac == 0,
        "AS_lowqual": ds.AS_lowqual,
        "AS_VQSR": hl.is_missing(ds.info["AS_VQSLOD"]),
    }
    snv_indel_expr = {'snv': hl.is_snp(ds.alleles[0], ds.alleles[1])}
    snv_indel_expr['indel'] = ~snv_indel_expr['snv']
    for var_type, score_cut in score_cutoffs.items():
        filters['AS_VQSR'] = filters['AS_VQSR'] | (snv_indel_expr[var_type] & (ds.info.AS_VQSLOD < score_cut.min_score))

    ds = ds.annotate(filters=add_filters_expr(filters=filters))

    # Drop unnecessary fields.
    logger.info('Dropping unnecessary fields...')
    summary_dict = {
        name.replace("expanded_", ""): ds[name]
        for name in [
            'colocated_variants',
            'subsets',
            'flags',
            'expanded_freq',
            'fafmax',
            'expanded_age_distribution',
            'vep',
            'filters',
            'quality_metrics',
        ]
    }
    ds = ds.annotate(variant_id=ds.variant_id, rsids=ds.rsids, summary=hl.struct(**summary_dict))
    ds = ds.select('variant_id', 'rsids', 'summary')
    ds = ds.rename({'summary': exome_or_genome})

    globals_dict = ds.globals.collect()[0]
    ds = ds.select_globals()
    ds = ds.annotate_globals(**{f"{exome_or_genome}_{k}": v for k, v in globals_dict.items()})

    return ds


def create_final_combined_faf_release(
    ht,
    contingency_table_ht: hl.Table,
    cmh_ht: hl.Table,
) -> hl.Table:
    """
    Create the final combined FAF release Table.

    :param ht: Table with joint exomes and genomes frequency and FAF information.
    :param contingency_table_ht: Table with contingency table test results to include
        on the final Table.
    :param cmh_ht: Table with Cochran–Mantel–Haenszel test results to include on the
        final Table.
    :return: Table with final combined FAF release information.
    """
    logger.info("Creating final combined FAF release Table...")

    def _get_struct_by_data_type(
        ht: hl.Table,
        annotation_expr=None,
        condition_func=None,
        postfix: str = "",
    ) -> hl.expr.StructExpression:
        """
        Group all annotations from the same data type into a struct.

        :param ht: Table with annotations to loop through.
        :param annotation_expr: StructExpression with annotations to loop through. If
            None, use ht.row.
        :param condition_func: Function that returns a BooleanExpression for a
            condition to check in addition to data type.
        :param postfix: String to add to the end of the new annotation following data
            type.
        :return: StructExpression of joint data type and the corresponding data type
            struct.
        """
        if condition_func is None:
            condition_func = lambda x: True
        if annotation_expr is None:
            annotation_expr = ht.row

        return hl.struct(
            **{
                f"{data_type}{postfix}": hl.struct(
                    **{
                        x.replace(f"{data_type}_", ""): ht[x]
                        for x in annotation_expr
                        if x.startswith(data_type) and condition_func(x)
                    },
                )
                for data_type in ['exomes', 'genomes', 'joint']
            },
        )

    # Group all the histograms into a single struct per data type.
    ht = ht.transmute(
        **_get_struct_by_data_type(ht, condition_func=lambda x: x.endswith("_hists"), postfix="_histograms"),
    )

    logger.info('Annotating capture and calling intervals.')
    capture_intervals_paths: list[str] = config_retrieve(['large_cohort', 'browser', 'capture_intervals'])
    capture_tables = [
        hl.import_locus_intervals(str(path), reference_genome=genome_build()) for path in capture_intervals_paths
    ]
    if len(capture_intervals_paths) > 1:
        # Deduplicate keys, keeping exactly one row for each unique key (key='interval').
        capture_intervals_ht = hl.Table.union(*capture_tables)
    else:
        capture_intervals_ht = capture_tables[0]

    calling_intervals_path = config_retrieve(['large_cohort', 'browser', 'calling_intervals'])
    calling_interval_padding = config_retrieve(['large_cohort', 'browser', 'calling_interval_padding'], default=150)
    calling_intervals_ht = hl.import_locus_intervals(str(calling_intervals_path), reference_genome=genome_build())

    intervals: list[tuple[str, hl.Table]] = [
        ('calling', adjust_interval_padding(calling_intervals_ht, calling_interval_padding)),
    ] + [
        ('capture', capture_intervals_ht),
    ]

    region_expr = {
        **{f"outside_{c}_region": hl.is_missing(i[ht.locus]) for c, i in intervals},
        **{
            f"not_called_in_{d}": hl.is_missing(ht[f"{d}_freq"]) | hl.is_missing(ht[f"{d}_freq"][1].AN)
            for d in ['exomes', 'genomes']
        },
    }

    # Get the stats expressions for the final Table.
    def _prep_stats_annotation(ht: hl.Table, stat_ht: hl.Table, stat_field: str):
        """
        Prepare stats expressions for the final Table.

        :param ht: Table to add stats to.
        :param stat_ht: Table with stats to add to the final Table.
        :param stat_field: Field to add to the final Table.
        :return: Dictionary of stats expressions for the final Table.
        """
        stat_expr = stat_ht[ht.key][stat_field]
        is_struct = False
        if isinstance(stat_expr, hl.expr.StructExpression):
            is_struct = True
            stat_expr = hl.array([stat_expr])

        stat_expr = stat_expr.map(
            lambda x: x.select(**{a: hl.or_missing(~hl.is_nan(x[a]), x[a]) for a in x if a != "gen_ancs"}),
        )

        if is_struct:
            stat_expr = stat_expr[0]

        return stat_expr

    # Prepare stats annotations for the final Table.
    ct_expr = _prep_stats_annotation(ht, contingency_table_ht, "contingency_table_test")
    cmh_expr = _prep_stats_annotation(ht, cmh_ht, "cochran_mantel_haenszel_test")
    anc_expr = cmh_ht[ht.key]["cochran_mantel_haenszel_test"].gen_ancs
    ct_union_expr = ct_expr[ht.joint_freq_index_dict[anc_expr[0] + "_adj"]]

    # Create a union stats annotation that uses the contingency table test if there is
    # only one genetic ancestry group, otherwise use the CMH test.
    union_expr = hl.if_else(
        hl.len(anc_expr) == 1,
        ct_union_expr.select("p_value").annotate(stat_test_name="contingency_table_test", gen_ancs=anc_expr),
        cmh_expr.select("p_value").annotate(stat_test_name="cochran_mantel_haenszel_test", gen_ancs=anc_expr),
    )
    stats_expr = {
        "contingency_table_test": ct_expr,
        "cochran_mantel_haenszel_test": cmh_expr,
        "stat_union": union_expr,
    }

    # Select the final columns for the Table, grouping all annotations from the same
    # data type into a struct.
    ht = ht.select(
        region_flags=hl.struct(**region_expr),
        **_get_struct_by_data_type(ht),
        freq_comparison_stats=hl.struct(**stats_expr),
    )

    # Group all global annotations from the same data type into a struct.
    ht = ht.transmute_globals(**_get_struct_by_data_type(ht, annotation_expr=ht.globals, postfix="_globals"))

    return ht


def prepare_gnomad_v4_variants_joint_frequency_helper(ds):
    globals = hl.eval(ds.globals)

    def joint_freq_index_key(subset=None, pop=None, sex=None, raw=False):
        parts = [s for s in [subset, pop, sex] if s is not None]
        parts.append("raw" if raw else "adj")
        return "_".join(parts)

    def freq_joint(ds, subset=None, pop=None, sex=None, raw=False):
        return ds.joint.freq[globals.joint_globals.freq_index_dict[joint_freq_index_key(subset, pop, sex, raw)]]

    flags = [
        hl.or_missing(ds.freq_comparison_stats.cochran_mantel_haenszel_test.p_value < 10e-4, "discrepant_frequencies"),
        hl.or_missing(ds.region_flags.not_called_in_exomes, "not_called_in_exomes"),
        hl.or_missing(ds.region_flags.not_called_in_genomes, "not_called_in_genomes"),
    ]

    ancestry_groups = set(m.get("pop", None) for m in globals.joint_globals.freq_meta)

    ds = ds.annotate(
        joint=hl.struct(
            freq=hl.struct(
                all=hl.struct(
                    ac=freq_joint(ds).AC,
                    ac_raw=freq_joint(ds, raw=True).AC,
                    an=freq_joint(ds).AN,
                    hemizygote_count=hl.if_else(
                        ds.locus.in_autosome_or_par(),
                        0,
                        hl.or_else(freq_joint(ds, sex="XY").AC, 0),
                    ),
                    homozygote_count=freq_joint(ds).homozygote_count,
                    ancestry_groups=[
                        hl.struct(
                            id="_".join(filter(bool, [pop, sex])),
                            ac=hl.or_else(freq_joint(ds, pop=pop, sex=sex).AC, 0),
                            an=hl.or_else(freq_joint(ds, pop=pop, sex=sex).AN, 0),
                            hemizygote_count=(
                                0
                                if sex == "XX"
                                else hl.if_else(
                                    ds.locus.in_autosome_or_par(),
                                    0,
                                    hl.or_else(freq_joint(ds, pop=pop, sex="XY").AC, 0),
                                )
                            ),
                            homozygote_count=hl.or_else(freq_joint(ds, pop=pop, sex=sex).homozygote_count, 0),
                        )
                        for pop, sex in list(itertools.product(ancestry_groups, [None, "XX", "XY"]))
                        + [(None, "XX"), (None, "XY")]
                    ],
                ),
            ),
            faf=ds.joint.faf,
            fafmax=ds.joint.fafmax,
            grpmax=ds.joint.grpmax,
            histograms=ds.joint.histograms,
            flags=hl.set(flags).filter(hl.is_defined),
            faf95_joint=hl.struct(
                grpmax=ds.joint.fafmax.faf95_max,
                grpmax_gen_anc=ds.joint.fafmax.faf95_max_gen_anc,
            ),
            faf99_joint=hl.struct(
                grpmax=ds.joint.fafmax.faf99_max,
                grpmax_gen_anc=ds.joint.fafmax.faf99_max_gen_anc,
            ),
            freq_comparison_stats=ds.freq_comparison_stats,
        ),
    )

    ds = ds.select(
        "joint",
    )

    return ds


def prepare_v4_variants(
    exome_ds_path: str,
    genome_ds_path: str,
    joint_freq_ht_path: str,
    contingency_ht_path: str,
    chm_ht_path: str,
    browser_outpath: str,
    exome_variants_outpath: str | None,
    genome_variants_outpath: str | None,
    joint_freq_outpath: str | None,
    transcript_table_path: str | None,
    mane_transcripts_path: str | None,
) -> hl.Table:

    # Process the score cutoffs.
    score_cutoffs = {
        seq_type: process_score_cutoffs(
            snv_bin_cutoff=config_retrieve(['large_cohort', 'browser', f'snp_bin_threshold_{seq_type}']),
            indel_bin_cutoff=config_retrieve(['large_cohort', 'browser', f'indel_bin_threshold_{seq_type}']),
            aggregated_bin_ht_path=config_retrieve(['large_cohort', 'browser', f'aggregated_bin_ht_path_{seq_type}']),
            snv_bin_id=config_retrieve(['large_cohort', 'browser', f'snp_bin_id_{seq_type}'], 'bin'),
            indel_bin_id=config_retrieve(['large_cohort', 'browser', f'indel_bin_id_{seq_type}'], 'bin'),
        )
        for seq_type in ['exome', 'genome']
    }

    exome_score_cutoffs = score_cutoffs['exome']
    genome_score_cutoffs = score_cutoffs['genome']

    if can_reuse(exome_variants_outpath):
        logger.info('Reusing existing exome_variants tables...')
        exome_variants = hl.read_table(exome_variants_outpath)
    else:
        logger.info('Preparing exome_variants tables...')
        exome_variants = prepare_gnomad_v4_variants_helper(exome_ds_path, 'exome', exome_score_cutoffs)
        exome_variants = exome_variants.checkpoint(exome_variants_outpath, overwrite=True)

    if can_reuse(genome_variants_outpath):
        logger.info('Reusing existing genome_variants tables...')
        genome_variants = hl.read_table(genome_variants_outpath)
    else:
        logger.info('Preparing genome_variants tables...')
        genome_variants = prepare_gnomad_v4_variants_helper(genome_ds_path, 'genome', genome_score_cutoffs)
        genome_variants = genome_variants.checkpoint(genome_variants_outpath, overwrite=True)

    # Add the key so that variant_id and rsid gets considered in the join
    genome_variants = genome_variants.key_by('locus', 'alleles', 'variant_id')
    exome_variants = exome_variants.key_by('locus', 'alleles', 'variant_id')

    # Drop rsids from exome_variants to avoid appending during join.
    exome_variants = exome_variants.drop('rsids')

    logger.info(f'Variants in exome: {exome_variants.count()}')
    logger.info(f'Variants in genome: {genome_variants.count()}')

    variants = exome_variants.join(genome_variants, how="outer")

    logger.info(f'Variants after join: {variants.count()}')

    shared_fields = [
        "vep",
    ]
    variants = variants.annotate(
        **{field: hl.or_else(variants.exome[field], variants.genome[field]) for field in shared_fields},
    )

    variants = variants.annotate(exome=variants.exome.drop(*shared_fields), genome=variants.genome.drop(*shared_fields))

    # Calculate colocated variants across exomes and genomes
    variants = variants.annotate(
        colocated_variants=hl.struct(
            all_exome=variants.exome.colocated_variants.all,
            all_genome=variants.genome.colocated_variants.all,
            **variants.exome.colocated_variants.drop("all"),
            **variants.genome.colocated_variants.drop("all"),
        ),
    )
    variants = variants.annotate(
        colocated_variants=variants.colocated_variants.annotate(
            all=hl.array(hl.set(hl.flatten([value for value in variants.colocated_variants.values()]))),
        ),
    )

    if can_reuse(joint_freq_outpath):
        joint_freq_ht = hl.read_table(joint_freq_outpath)
    else:
        logger.info('Preparing joint frequency table...')
        joint_freq_ht = hl.read_table(joint_freq_ht_path)
        contingency_table_ht = hl.read_table(contingency_ht_path)
        chm_ht = hl.read_table(chm_ht_path)

        joint_freq_ht = create_final_combined_faf_release(joint_freq_ht, contingency_table_ht, chm_ht)
        joint_freq_ht = joint_freq_ht.checkpoint(joint_freq_outpath, overwrite=True)

    joint_freq_ht = prepare_gnomad_v4_variants_joint_frequency_helper(joint_freq_ht)

    logger.info('Annotating variants with joint frequencies...')
    variants = variants.annotate(**joint_freq_ht[variants.locus, variants.alleles])

    logger.info('Annotating transcript consequences...')
    base_transcripts_grch38 = hl.read_table(transcript_table_path)

    transcripts = extract_transcripts(base_transcripts_grch38)

    variants = annotate_transcript_consequences(
        variants,
        transcripts,
        mane_transcripts_path,
    )

    logger.info('Re-keying variants by locus and alleles...')
    variants = variants.key_by('locus', 'alleles')

    # Return the final browser table.
    logger.info('Writing the final browser table...')
    variants.write(browser_outpath, overwrite=True)

    return variants
