import itertools
import logging

import hail as hl

from gnomad.utils.filtering import add_filters_expr

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
                    `5utr_consequence`: str,
                    `5utr_annotation`: dict<str, struct {
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


def annotate_transcript_consequences(ds, transcripts_path, mane_transcripts_path=None):

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

    transcripts = hl.read_table(transcripts_path)

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


def prepare_gnomad_v4_variants_helper(ds_path: str | None, exome_or_genome: str) -> hl.Table:

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

    logging.info('Preparing the base frequencies table...')
    ds = hl.read_table(ds_path)

    subsets = set(m.get("subset", None) for m in ds.freq_meta.collect()[0])

    ds = ds.annotate(variant_id=_variant_id(ds.locus, ds.alleles))
    ds = ds.rename({"rsid": "rsids"})

    # Identify co-located variants.
    logging.info('Identifying co-located variants...')
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
    logging.info('Expanding the frequencies struct...')

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
    logging.info('Identifying the subsets that contain each variant...')
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
    logging.info('Summarising the age distribution for variant carriers...')
    # ds = ds.annotate(expanded_age_distribution=hl.struct(het=ds.age_hist_het, hom=ds.age_hist_hom))

    # Summarise the site quality metrics.
    logging.info('Summarising the site quality metrics...')
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
    logging.info('Restructuring the region flags...')
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
    logging.info('Adding the variant filter flags...')

    inbreeding_coeff_cutoff = -0.3
    filters = {
        "InbreedingCoeff": ds.inbreeding_coeff[0] < inbreeding_coeff_cutoff,
        "AC0": ds.expanded_freq.all.ac == 0,
        "AS_lowqual": ds.AS_lowqual,
        "AS_VQSR": ds.as_vqsr_filters != "PASS",
    }
    ds = ds.annotate(filters=add_filters_expr(filters=filters))

    # Drop unnecessary fields.
    logging.info('Dropping unnecessary fields...')
    summary_dict = {
        name.replace("expanded_", ""): ds[name]
        for name in [
            'colocated_variants',
            'subsets',
            'flags',
            'expanded_freq',
            'fafmax',
            # 'expanded_age_distribution',
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


def prepare_v4_variants(
    exome_ds_path: str,
    genome_ds_path: str,
    browser_outpath: str,
    exome_variants_outpath: str | None,
    genome_variants_outpath: str | None,
    transcript_table_path: str | None,
    mane_transcripts_path: str | None,
) -> hl.Table:

    # Generate the browser output tables for each data type.
    exome_variants = prepare_gnomad_v4_variants_helper(exome_ds_path, 'exome')
    genome_variants = prepare_gnomad_v4_variants_helper(genome_ds_path, 'genome')

    # checkpoint datasets
    exome_variants = exome_variants.checkpoint(exome_variants_outpath, overwrite=True)
    genome_variants = genome_variants.checkpoint(genome_variants_outpath, overwrite=True)

    # Add the key so that variant_id and rsid gets considered in the join
    genome_variants = genome_variants.key_by('locus', 'alleles', 'variant_id', 'rsids')
    exome_variants = exome_variants.key_by('locus', 'alleles', 'variant_id', 'rsids')
    variants = exome_variants.join(genome_variants, how="outer")

    shared_fields = [
        # "lcr",
        # "nonpar",
        # "rsids",
        # "segdup",
        "vep",
        # "in_silico_predictors",
        # "variant_id",
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

    variants = variants.checkpoint(browser_outpath, overwrite=True)

    logging.info('Annotating transcript consequences...')
    variants = annotate_transcript_consequences(
        variants,
        transcript_table_path,
        mane_transcripts_path,
    )

    # Return the final browser table.
    variants.write(browser_outpath, overwrite=True)

    return variants
