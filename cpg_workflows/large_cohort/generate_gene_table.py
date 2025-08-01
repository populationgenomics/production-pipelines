"""
This script implements a pipeline for generating a public gene annotation table for GRCh38, suitable for browser and downstream analyses.

It is based on the gnomAD browser's gene pipeline:
https://github.com/broadinstitute/gnomad-browser/blob/main/data-pipeline/src/data_pipeline/pipelines/genes.py
with helper functions adapted from:
https://github.com/broadinstitute/gnomad-browser/blob/main/data-pipeline/src/data_pipeline/data_types/gene.py

Main steps:
- Imports and processes GENCODE gene, transcript, and exon annotations, and HGNC gene labels.
- Extracts canonical and MANE Select transcripts from variant frequency tables and MANE summary files.
- Annotates genes and transcripts with preferred transcript IDs, RefSeq IDs, and other metadata.
- Merges overlapping exons and organizes exon features for browser display.
- Filters out pseudoautosomal region Y (PAR_Y) genes.
- Checkpoints intermediate results for reproducibility and downstream use.
- Prepares the final gene table for public release, including global annotations.
"""

import logging

import pandas as pd

import hail as hl

from cpg_utils.config import reference_path
from cpg_workflows.utils import can_reuse

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

HGNC_COLUMNS = [
    "gd_hgnc_id",
    "gd_app_sym",
    "gd_app_name",
    "gd_prev_sym",
    "gd_aliases",
    "gd_pub_eg_id",
    "gd_pub_ensembl_id",
    "md_eg_id",
    "md_ensembl_id",
    "md_mim_id",
]


def import_mane_select_transcripts(path, version: str = 'v1.4') -> hl.Table:

    ds = hl.import_table(path, force=True)

    ds = ds.annotate_globals(version=version)

    ds = ds.filter(ds.MANE_status == "MANE Select")

    ds = ds.select(
        gene_id=ds.Ensembl_Gene.split("\\.")[0],
        matched_gene_version=ds.Ensembl_Gene.split("\\.")[1],
        ensembl_id=ds.Ensembl_nuc.split("\\.")[0],
        ensembl_version=ds.Ensembl_nuc.split("\\.")[1],
        refseq_id=ds.RefSeq_nuc.split("\\.")[0],
        refseq_version=ds.RefSeq_nuc.split("\\.")[1],
    )

    return ds.key_by("gene_id")


def get_canonical_transcripts(**sites_tables) -> hl.Table:
    canonical_transcripts: set[tuple[str, str]] = set()
    for sites_table in sites_tables.values():
        table_canonical_transcripts = sites_table.aggregate(
            hl.agg.explode(
                lambda csq: hl.agg.collect_as_set((csq.gene_id, csq.transcript_id)),
                sites_table.vep.transcript_consequences.filter(lambda csq: csq.canonical == 1),
            ),
        )
        canonical_transcripts = canonical_transcripts.union(table_canonical_transcripts)

    canonical_transcripts_ht: hl.Table = hl.Table.from_pandas(
        pd.DataFrame(
            {"gene_id": gene_id, "canonical_transcript_id": canonical_transcript_id}
            for gene_id, canonical_transcript_id in canonical_transcripts
        ),
        key="gene_id",
    )

    return canonical_transcripts_ht.repartition(32, shuffle=True)


def normalized_contig(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.rbind(hl.str(contig).replace("^chr", ""), lambda c: hl.if_else(c == "MT", "M", c))


def contig_number(contig: hl.expr.StringExpression) -> hl.expr.Int32Expression:
    return hl.bind(
        lambda contig: (
            hl.case().when(contig == "X", 23).when(contig == "Y", 24).when(contig == "M", 25).default(hl.int(contig))
        ),
        normalized_contig(contig),
    )


def x_position(locus: hl.expr.LocusExpression) -> hl.expr.Int64Expression:
    return hl.int64(contig_number(locus.contig)) * 1_000_000_000 + locus.position


def prepare_genes(gencode_path, hgnc_path, reference_genome):
    genes = import_gencode(gencode_path, reference_genome)

    hgnc = import_hgnc(hgnc_path)
    hgnc = hgnc.filter(hl.is_defined(hgnc.ensembl_id)).key_by("ensembl_id")
    genes = genes.annotate(**hgnc[genes.gene_id])
    # If a symbol was not present in HGNC data, use the symbol from GENCODE
    genes = genes.annotate(symbol=hl.or_else(genes.symbol, genes.gencode_symbol))

    genes = genes.annotate(
        symbol_upper_case=genes.symbol.upper(),
        search_terms=hl.empty_set(hl.tstr)
        .add(genes.symbol)
        .add(genes.gencode_symbol)
        .union(genes.previous_symbols)
        .union(genes.alias_symbols)
        .map(lambda s: s.upper()),
    )

    genes = genes.annotate(
        reference_genome=reference_genome,
        transcripts=genes.transcripts.map(lambda transcript: transcript.annotate(reference_genome=reference_genome)),
    )

    chip_genes = {"ENSG00000171456", "ENSG00000119772", "ENSG00000168769"}
    genes = genes.annotate(
        flags=hl.set([hl.or_missing(hl.set(chip_genes).contains(genes.gene_id), "chip")]).filter(hl.is_defined),
    )

    return genes


def merge_overlapping_exons(regions):
    return hl.if_else(
        hl.len(regions) > 1,
        hl.rbind(
            hl.sorted(regions, lambda region: region.start),
            lambda sorted_regions: sorted_regions[1:].fold(
                lambda acc, region: hl.if_else(
                    region.start <= acc[-1].stop + 1,
                    acc[:-1].append(
                        acc[-1].annotate(
                            stop=hl.max(region.stop, acc[-1].stop),
                            xstop=hl.max(region.xstop, acc[-1].xstop),
                        ),
                    ),
                    acc.append(region),
                ),
                [sorted_regions[0]],
            ),
        ),
        regions,
    )


###############################################
# Exons                                       #
###############################################


def get_exons(gencode):
    """
    Filter GENCODE table to exons and format fields.
    """
    exons = gencode.filter(hl.set(["exon", "CDS", "UTR"]).contains(gencode.feature))
    exons = exons.select(
        feature_type=exons.feature,
        transcript_id=exons.transcript_id.split("\\.")[0],
        transcript_version=exons.transcript_id.split("\\.")[1],
        gene_id=exons.gene_id.split("\\.")[0],
        gene_version=exons.gene_id.split("\\.")[1],
        chrom=normalized_contig(exons.interval.start.contig),
        strand=exons.strand,
        start=exons.interval.start.position,
        stop=exons.interval.end.position,
        xstart=x_position(exons.interval.start),
        xstop=x_position(exons.interval.end),
    )

    return exons


###############################################
# Genes                                       #
###############################################


def get_genes(gencode):
    """
    Filter GENCODE table to genes and format fields.
    """
    genes = gencode.filter(gencode.feature == "gene")
    genes = genes.select(
        gene_id=genes.gene_id.split("\\.")[0],
        gene_version=genes.gene_id.split("\\.")[1],
        gencode_symbol=genes.gene_name,
        chrom=normalized_contig(genes.interval.start.contig),
        strand=genes.strand,
        start=genes.interval.start.position,
        stop=genes.interval.end.position,
        xstart=x_position(genes.interval.start),
        xstop=x_position(genes.interval.end),
    )

    genes = genes.annotate()

    genes = genes.key_by(genes.gene_id)

    return genes


def collect_gene_exons(gene_exons):
    # There are 3 feature types in the exons collection: "CDS", "UTR", and "exon".
    # There are "exon" regions that cover the "CDS" and "UTR" regions and also
    # some (non-coding) transcripts that contain only "exon" regions.
    # This filters the "exon" regions to only those that are in non-coding transcripts.
    #
    # This makes the UI for selecting visible regions easier, since it can filter
    # on "CDS" or "UTR" feature type without having to also filter out the "exon" regions
    # that duplicate the "CDS" and "UTR" regions.

    non_coding_transcript_exons = hl.bind(
        lambda coding_transcripts: gene_exons.filter(lambda exon: ~coding_transcripts.contains(exon.transcript_id)),
        hl.set(
            gene_exons.filter(lambda exon: (exon.feature_type == "CDS") | (exon.feature_type == "UTR")).map(
                lambda exon: exon.transcript_id,
            ),
        ),
    )

    exons = (
        merge_overlapping_exons(gene_exons.filter(lambda exon: exon.feature_type == "CDS"))
        .extend(merge_overlapping_exons(gene_exons.filter(lambda exon: exon.feature_type == "UTR")))
        .extend(merge_overlapping_exons(non_coding_transcript_exons))
    )

    exons = exons.map(lambda exon: exon.select("feature_type", "start", "stop", "xstart", "xstop"))

    return exons


def reject_par_y_genes(genes: hl.Table) -> hl.Table:
    genes = genes.filter(genes.gene_version.endswith("_PAR_Y") == hl.literal(False))
    return genes


###############################################
# Transcripts                                 #
###############################################


def get_transcripts(gencode):
    """
    Filter GENCODE table to transcripts and format fields.
    """
    transcripts = gencode.filter(gencode.feature == "transcript")
    transcripts = transcripts.select(
        transcript_id=transcripts.transcript_id.split("\\.")[0],
        transcript_version=transcripts.transcript_id.split("\\.")[1],
        gene_id=transcripts.gene_id.split("\\.")[0],
        gene_version=transcripts.gene_id.split("\\.")[1],
        chrom=normalized_contig(transcripts.interval.start.contig),
        strand=transcripts.strand,
        start=transcripts.interval.start.position,
        stop=transcripts.interval.end.position,
        xstart=x_position(transcripts.interval.start),
        xstop=x_position(transcripts.interval.end),
    )

    return transcripts.key_by(transcripts.transcript_id)


def collect_transcript_exons(transcript_exons):
    # There are 3 feature types in the exons collection: "CDS", "UTR", and "exon".
    # There are "exon" regions that cover the "CDS" and "UTR" regions and also
    # some (non-coding) transcripts that contain only "exon" regions.
    # This filters the "exon" regions to only those that are in non-coding transcripts.
    #
    # This makes the UI for selecting visible regions easier, since it can filter
    # on "CDS" or "UTR" feature type without having to also filter out the "exon" regions
    # that duplicate the "CDS" and "UTR" regions.

    is_coding = transcript_exons.any(lambda exon: (exon.feature_type == "CDS") | (exon.feature_type == "UTR"))

    exons = hl.if_else(is_coding, transcript_exons.filter(lambda exon: exon.feature_type != "exon"), transcript_exons)

    exons = exons.map(lambda exon: exon.select("feature_type", "start", "stop", "xstart", "xstop"))

    return exons


###############################################
# Main                                        #
###############################################


def import_gencode(path, reference_genome):
    gencode = hl.experimental.import_gtf(path, force=True, reference_genome=reference_genome)
    gencode = gencode.repartition(2000, shuffle=True)
    gencode = gencode.cache()

    # Extract genes and transcripts
    genes = get_genes(gencode)
    transcripts = get_transcripts(gencode)

    # Annotate genes/transcripts with their exons
    exons = get_exons(gencode)
    exons = exons.cache()

    gene_exons = exons.group_by(exons.gene_id, exons.gene_version).aggregate(exons=hl.agg.collect(exons.row_value))
    genes = genes.annotate(exons=collect_gene_exons(gene_exons[genes.gene_id, genes.gene_version].exons))

    transcript_exons = exons.group_by(exons.transcript_id, exons.transcript_version).aggregate(
        exons=hl.agg.collect(exons.row_value),
    )
    transcripts = transcripts.annotate(
        exons=collect_transcript_exons(
            transcript_exons[transcripts.transcript_id, transcripts.transcript_version].exons,
        ),
    )

    # Annotate genes with their transcripts
    gene_transcripts = transcripts.key_by()
    gene_transcripts = gene_transcripts.group_by(gene_transcripts.gene_id, gene_transcripts.gene_version).aggregate(
        transcripts=hl.agg.collect(gene_transcripts.row_value),
    )
    genes = genes.annotate(transcripts=gene_transcripts[genes.gene_id, genes.gene_version].transcripts)

    return genes


def import_hgnc(path):
    ds = hl.import_table(path, missing="")

    ds = ds.select(
        hgnc_id=ds["HGNC ID"],
        symbol=ds["Approved symbol"],
        name=ds["Approved name"],
        previous_symbols=ds["Previous symbols"],
        alias_symbols=ds["Alias symbols"],
        omim_id=ds["OMIM ID(supplied by OMIM)"],
        ensembl_id=hl.or_else(ds["Ensembl gene ID"], ds["Ensembl ID(supplied by Ensembl)"]),
        ncbi_id=hl.or_else(ds["NCBI Gene ID"], ds["NCBI Gene ID(supplied by NCBI)"]),
    )

    ds = ds.annotate(
        previous_symbols=hl.set(ds.previous_symbols.split(",").map(lambda s: s.strip())),
        alias_symbols=hl.set(ds.alias_symbols.split(",").map(lambda s: s.strip())),
    )

    return ds


def prepare_gene_table_for_release(ds: hl.Table, keep_mane_version_global_annotation: bool = False) -> hl.Table:
    if keep_mane_version_global_annotation:
        globals_dict = ds.index_globals()
        ds = ds.select_globals(mane_select_version=globals_dict["annotations"]["mane_select_transcript"]["version"])
    else:
        ds = ds.select_globals()

    ds = ds.repartition(50)
    return ds


def annotate_table(ds: hl.Table, join_on=None, **annotation_tables) -> hl.Table:

    for annotation_key, annotation_table in annotation_tables.items():

        ds = ds.annotate_globals(
            annotations=getattr(ds.globals, "annotations", hl.struct()).annotate(
                **{annotation_key: hl.eval(annotation_table.globals)},
            ),
        )

        # If the table has only one non-key field, add that field's value directly to the table.
        # Otherwise, add a struct field containing the annotation table's value.
        if len(annotation_table.row_value) == 1:
            ds = ds.annotate(**annotation_table[ds[join_on] if join_on else ds.key])
        else:
            ds = ds.annotate(**{annotation_key: annotation_table[ds[join_on] if join_on else ds.key]})

    return ds


def annotate_gene_transcripts_with_refseq_id(genes: hl.Table, mane_select_transcripts: hl.Table) -> hl.Table:

    ensembl_to_refseq_map = {}
    for transcript in mane_select_transcripts.collect():
        ensembl_to_refseq_map[transcript.ensembl_id] = {
            transcript.ensembl_version: hl.Struct(
                refseq_id=transcript.refseq_id,
                refseq_version=transcript.refseq_version,
            ),
        }

    ensembl_to_refseq_map = hl.literal(ensembl_to_refseq_map)

    genes = genes.annotate(
        transcripts=genes.transcripts.map(
            lambda transcript: transcript.annotate(
                **ensembl_to_refseq_map.get(
                    transcript.transcript_id,
                    hl.empty_dict(hl.tstr, hl.tstruct(refseq_id=hl.tstr, refseq_version=hl.tstr)),
                ).get(
                    transcript.transcript_version,
                    hl.struct(refseq_id=hl.null(hl.tstr), refseq_version=hl.null(hl.tstr)),
                ),
            ),
        ),
    )

    return genes


def annotate_with_preferred_transcript(ds: hl.Table) -> hl.Table:

    if "mane_select_transcript" in ds.row:
        preferred_transcript_id = hl.or_else(ds.mane_select_transcript.ensembl_id, ds.canonical_transcript_id)
        preferred_transcript_source = "mane_select"
    else:
        preferred_transcript_id = ds.canonical_transcript_id
        preferred_transcript_source = "ensembl_canonical"

    return ds.annotate(
        preferred_transcript_id=preferred_transcript_id,
        preferred_transcript_source=preferred_transcript_source,
    )


def run(
    genome_freq_ht_path: str | None,
    exome_freq_ht_path: str | None,
    tmp_prefix: str,
    step_six_output_path: str,
    mane_select_transcripts_ht_path: str,
    output_path: str,
):
    if can_reuse(step_six_output_path):
        logger.info(f"Reusing step 6 output from {step_six_output_path}")
        annotate_grch38_genes_step_6 = hl.read_table(str(step_six_output_path))
    else:
        freq_ht = hl.read_table(str(genome_freq_ht_path))
        exome_freq_ht = hl.read_table(str(exome_freq_ht_path))

        logger.info("Importing MANE Select transcripts...")
        mane_select_transcripts = import_mane_select_transcripts(reference_path('mane_1.4/summary'))
        mane_select_transcripts = mane_select_transcripts.checkpoint(
            mane_select_transcripts_ht_path,
            overwrite=True,
        )

        logger.info("Extracting canonical transcripts from frequency table...")
        canonical_transcripts_grch38 = get_canonical_transcripts(genomes=freq_ht, exomes=exome_freq_ht)

        logger.info("Preparing base gene table from GENCODE and HGNC...")
        genes_grch38_base: hl.Table = prepare_genes(
            gencode_path=reference_path('gencode_v44'),
            hgnc_path=reference_path('hgnc_labels'),
            reference_genome='GRCh38',
        )

        genes_grch38_base = genes_grch38_base.checkpoint(tmp_prefix + '/genes_grch38_base.ht', overwrite=True)

        logger.info("Annotating gene table with canonical and MANE Select transcripts...")
        annotate_grch38_genes_step_1 = annotate_table(
            genes_grch38_base,
            canonical_transcript=canonical_transcripts_grch38,
            mane_select_transcript=mane_select_transcripts,
        )

        annotate_grch38_genes_step_1 = annotate_grch38_genes_step_1.checkpoint(
            tmp_prefix + '/genes_grch38_annotated_step_1.ht',
            overwrite=True,
        )

        # SKIPPING STEP 2: ANNOTATION WITH GTEX

        logger.info("Annotating gene transcripts with RefSeq IDs...")
        genes_grch38_annotated_3 = annotate_gene_transcripts_with_refseq_id(
            annotate_grch38_genes_step_1,
            mane_select_transcripts,
        )

        genes_grch38_annotated_3 = genes_grch38_annotated_3.checkpoint(
            tmp_prefix + '/genes_grch38_annotated_step_3.ht',
            overwrite=True,
        )

        logger.info("Annotating with preferred transcript...")
        annotate_grch38_genes_step_4 = annotate_with_preferred_transcript(
            genes_grch38_annotated_3,
        )

        annotate_grch38_genes_step_4 = annotate_grch38_genes_step_4.checkpoint(
            tmp_prefix + '/genes_grch38_annotated_step_4.ht',
            overwrite=True,
        )

        # SKIPPING STEP 5: ANNOTATING WITH CONSTRAINT

        logger.info("Filtering out PAR_Y genes...")
        annotate_grch38_genes_step_6 = reject_par_y_genes(
            annotate_grch38_genes_step_4,
        )

        logger.info("Checkpointing step 6 output for PrepareBrowserTable...")
        # Step 6 output is used in PrepareBrowserTable to extract transcripts for browser so is written to non-tmp bucket
        annotate_grch38_genes_step_6 = annotate_grch38_genes_step_6.checkpoint(
            step_six_output_path,
            overwrite=True,
        )

    # SKIPPING REMOVING CONSTRAINT ANNOTATION BECAUSE WE NEVER ANNOTATED WITH IT

    logger.info("Preparing gene table for public release...")
    prepare_grch38_genes_table_for_public_release = prepare_gene_table_for_release(
        annotate_grch38_genes_step_6,
        True,
    )

    prepare_grch38_genes_table_for_public_release = prepare_grch38_genes_table_for_public_release.checkpoint(
        output_path,
        overwrite=True,
    )
