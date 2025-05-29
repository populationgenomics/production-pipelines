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


def prepare_gnomad_v4_variants_helper(ds_path: str | None, exome_or_genome: str) -> hl.Table:

    if ds_path is None:
        ds = hl.Table.parallelize(hl.literal([], "".join(EMPTY_TABLE_CONFIG.split()))).key_by('locus', 'alleles')
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
    ds = ds.annotate(expanded_age_distribution=hl.struct(het=ds.age_hist_het, hom=ds.age_hist_hom))

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
            'expanded_age_distribution',
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

    # Return the final browser table.
    variants.write(browser_outpath, overwrite=True)

    return variants
