"""
Hail Query functions for seqr loader.
"""

import logging

import hail as hl

from cpg_utils.config import config_retrieve, get_config, reference_path
from cpg_utils.hail_batch import genome_build
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.large_cohort.load_vqsr import load_vqsr
from cpg_workflows.utils import checkpoint_hail
from hail_scripts.computed_fields import variant_id, vep


def annotate_cohort(
    vcf_path: str,
    out_mt_path: str,
    vep_ht_path: str,
    site_only_vqsr_vcf_path=None,
    checkpoint_prefix=None,
    long_read=False,
):
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR
    annotations.
    """

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    # we satisfy this requirement by calling `init_batch()` in the query_command wrapper
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    # tune the logger correctly
    logging.getLogger().setLevel(logging.INFO)

    # hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/permissions.20issues/near/398711114
    # don't override the block size, as it explodes the number of partitions when processing TB+ datasets
    # Each partition comes with some computational overhead, it's to be seen whether the standard block size
    # is viable for QOB in large datasets... Requires a test
    mt = hl.import_vcf(
        vcf_path,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
        force_bgz=True,
        array_elements_required=False,
    )
    logging.info(f'Imported VCF {vcf_path} as {mt.n_partitions()} partitions')

    # Annotate VEP. Do it before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(vep_ht_path)
    logging.info(
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

    if site_only_vqsr_vcf_path:
        vqsr_ht = load_vqsr(site_only_vqsr_vcf_path)
        vqsr_ht = checkpoint_hail(vqsr_ht, 'vqsr.ht', checkpoint_prefix)

        logging.info('Adding VQSR annotations into the Matrix Table')
        mt = mt.annotate_globals(**vqsr_ht.index_globals())
        mt = mt.annotate_rows(
            # vqsr_ht has info annotation split by allele, plus the new AS-VQSR annotations
            info=vqsr_ht[mt.row_key].info,
            filters=mt.filters.union(vqsr_ht[mt.row_key].filters).filter(lambda val: val != 'PASS'),
        )
        mt = checkpoint_hail(mt, 'mt-vep-split-vqsr.mt', checkpoint_prefix)

    ref_ht = hl.read_table(reference_path('seqr_combined_reference_data'))
    clinvar_ht = hl.read_table(reference_path('seqr_clinvar'))

    logging.info('Annotating with seqr-loader fields: round 1')

    # split the AC/AF attributes into separate entries, overwriting the array in INFO
    # these elements become a 1-element array
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AF=[mt.info.AF[mt.a_index - 1]],
            AC=[mt.info.AC[mt.a_index - 1]],
        ),
    )

    logging.info('Annotating with clinvar and munging annotation fields')
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
    logging.info('Adding GRCh37 coords')
    liftover_path = reference_path('liftover_38_to_37')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_path, rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))

    # only remove InbreedingCoeff if present (post-VQSR)
    if 'InbreedingCoeff' in mt.info:
        mt = mt.annotate_rows(info=mt.info.drop('InbreedingCoeff'))

    logging.info(
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
        sourceFilePath=vcf_path,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
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

    logging.info('Done:')
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Written final matrix table into {out_mt_path}')


def subset_mt_to_samples(mt_path: str, sample_ids: list[str], out_mt_path: str, exclusion_file: str | None = None):
    """
    Subset the MatrixTable to the provided list of samples and to variants present
    in those samples

    Args:
        mt_path (str): cohort-level matrix table from VCF.
        sample_ids (list[str]): samples to take from the matrix table.
        out_mt_path (str): path to write the result.
        exclusion_file (str, optional): path to a file containing samples to remove from the
                                        subset prior to extracting
    """

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    # we satisfy this requirement by calling `init_batch()` in the query_command wrapper
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    logging.basicConfig(level=logging.INFO)

    mt = hl.read_matrix_table(mt_path)

    unique_sids: set[str] = set(sample_ids)

    # if an exclusion file was passed, remove the samples from the subset
    # this executes in a query command, by execution time the file should exist
    if exclusion_file:
        with hl.hadoop_open(exclusion_file) as f:
            exclusion_ids = set(f.read().splitlines())
        logging.info(f'Excluding {len(exclusion_ids)} samples from the subset')
        unique_sids -= exclusion_ids

    mt_sample_ids = set(mt.s.collect())

    if sample_ids_not_in_mt := unique_sids - mt_sample_ids:
        raise ValueError(
            f'Found {len(sample_ids_not_in_mt)}/{len(unique_sids)} IDs in the requested subset not in the callset.\n'
            f'IDs that aren\'t in the callset: {sample_ids_not_in_mt}\n'
            f'All callset sample IDs: {mt_sample_ids}',
        )

    logging.info(f'Found {len(mt_sample_ids)} samples in mt, subsetting to {len(unique_sids)} samples.')

    mt = mt.filter_cols(hl.literal(unique_sids).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Finished subsetting to {len(unique_sids)} samples, written to {out_mt_path}.')


def vcf_from_mt_subset(mt_path: str, out_vcf_path: str):
    """
    Read the MT in, and write out to a VCF
    If we wanted to translate sample IDs to external samples
    then we could do that here, otherwise rely on VCF re-heading

    Args:
        mt_path (str): path of the single-dataset MT to read in
        out_vcf_path (str): path of the vcf.bgz to generate
    """

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    # we satisfy this requirement by calling `init_batch()` in the query_command wrapper
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    mt = hl.read_matrix_table(mt_path)
    logging.info(f'Dataset MT dimensions: {mt.count()}')
    hl.export_vcf(mt, out_vcf_path, tabix=True)
    logging.info(f'Written {out_vcf_path}')


def annotate_dataset_mt(mt_path: str, out_mt_path: str):
    """
    Add dataset-level annotations.
    """

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    # we satisfy this requirement by calling `init_batch()` in the query_command wrapper
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
        override_jar_spec(jar_spec)

    mt = hl.read_matrix_table(mt_path)

    # Convert the mt genotype entries into num_alt, gq, ab, dp, and sample_id.
    is_called = hl.is_defined(mt.GT)
    genotype_fields = {
        'num_alt': hl.if_else(is_called, mt.GT.n_alt_alleles(), -1),
        'gq': hl.if_else(is_called, mt.GQ, hl.null(hl.tint)),
        'ab': hl.bind(
            lambda total: hl.if_else(
                is_called & (total != 0) & (hl.len(mt.AD) > 1),
                hl.float(mt.AD[1] / total),
                hl.missing(hl.tfloat),
            ),
            hl.sum(mt.AD),
        ),
        'dp': hl.if_else(is_called, hl.int(hl.min(mt.DP, 32000)), hl.missing(hl.tfloat)),
        'sample_id': mt.s,
    }
    logging.info('Annotating genotypes')
    mt = mt.annotate_rows(genotypes=hl.agg.collect(hl.struct(**genotype_fields)))

    def _genotype_filter_samples(fn):
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    # 2022-07-28 mfranklin: Initially the code looked like:
    #           {**_genotype_filter_samples(lambda g: g.num_alt == i) for i in ...}
    #   except the lambda definition doesn't bind the loop variable i in this scope
    #   instead let's define the filters as functions, and wrap them in a decorator
    #   that captures the value of i.

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
        return (g.gq >= i) & (g.gq < i + 5)

    @_capture_i_decorator
    def _filter_samples_ab(i, g):
        return (g.num_alt == 1) & ((g.ab * 100) >= i) & ((g.ab * 100) < i + 5)

    mt = mt.annotate_rows(
        samples_no_call=_genotype_filter_samples(lambda g: g.num_alt == -1),
        samples_num_alt=hl.struct(**{('%i' % i): _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
        samples_gq=hl.struct(
            **{('%i_to_%i' % (i, i + 5)): _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 95, 5)},
        ),
        samples_ab=hl.struct(
            **{'%i_to_%i' % (i, i + 5): _genotype_filter_samples(_filter_samples_ab(i)) for i in range(0, 45, 5)},
        ),
    )
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Written {out_mt_path}')
