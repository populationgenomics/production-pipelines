"""
Hail Query functions for seqr loader.
"""

import logging
import os

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, genome_build
from hail_scripts.computed_fields import vep, variant_id

from cpg_workflows.large_cohort.load_vqsr import load_vqsr
from cpg_workflows.utils import can_reuse


def annotate_cohort(
    vcf_path,
    out_mt_path,
    vep_ht_path,
    site_only_vqsr_vcf_path=None,
    checkpoint_prefix=None,
):
    """
    Convert VCF to matrix table, annotate for Seqr Loader, add VEP and VQSR
    annotations.
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
            if can_reuse(path):
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
    logging.info(f'Importing VCF {vcf_path}')

    logging.info(f'Loading VEP Table from {vep_ht_path}')
    # Annotate VEP. Do ti before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(str(vep_ht_path))
    logging.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # Splitting multi-allelics. We do not handle AS info fields here - we handle
    # them when loading VQSR instead, and populate entire "info" from VQSR.
    mt = hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )
    mt = _checkpoint(mt, 'mt-vep-split.mt')

    if site_only_vqsr_vcf_path:
        vqsr_ht = load_vqsr(site_only_vqsr_vcf_path)
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

    ref_ht = hl.read_table(str(reference_path('seqr_combined_reference_data')))
    clinvar_ht = hl.read_table(str(reference_path('seqr_clinvar')))

    logging.info('Annotating with seqr-loader fields: round 1')

    # don't fail if the AC/AF attributes are an inappropriate type
    for attr in ['AC', 'AF']:
        if not isinstance(mt.info[attr], hl.ArrayExpression):
            mt = mt.annotate_rows(info=mt.info.annotate(**{attr: [mt.info[attr]]}))

    mt = mt.annotate_rows(
        AC=mt.info.AC[mt.a_index - 1],
        AF=mt.info.AF[mt.a_index - 1],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=variant_id.get_expr_for_variant_ids(
            mt.locus_old, mt.alleles_old
        ),
        sortedTranscriptConsequences=vep.get_expr_for_vep_sorted_transcript_consequences_array(
            mt.vep
        ),
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
    mt = _checkpoint(mt, 'mt-vep-split-vqsr-round1.mt')

    logging.info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)'
    )
    mt = mt.annotate_rows(
        domains=vep.get_expr_for_vep_protein_domains_set_from_sorted(
            mt.sortedTranscriptConsequences
        ),
        transcriptConsequenceTerms=vep.get_expr_for_vep_consequence_terms_set(
            mt.sortedTranscriptConsequences
        ),
        transcriptIds=vep.get_expr_for_vep_transcript_ids_set(
            mt.sortedTranscriptConsequences
        ),
        mainTranscript=vep.get_expr_for_worst_transcript_consequence_annotations_struct(
            mt.sortedTranscriptConsequences
        ),
        geneIds=vep.get_expr_for_vep_gene_ids_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep.get_expr_for_vep_gene_ids_set(
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
    mt.write(str(out_mt_path), overwrite=True)
    logging.info(f'Written final matrix table into {out_mt_path}')


def subset_mt_to_samples(mt_path, sample_ids, out_mt_path):
    """
    Subset the MatrixTable to the provided list of samples and to variants present
    in those samples
    @param mt_path: cohort-level matrix table from VCF.
    @param sample_ids: samples to take from the matrix table.
    @param out_mt_path: path to write the result.
    """
    mt = hl.read_matrix_table(str(mt_path))

    sample_ids = set(sample_ids)
    mt_sample_ids = set(mt.s.collect())

    sample_ids_not_in_mt = sample_ids - mt_sample_ids
    if sample_ids_not_in_mt:
        raise Exception(
            f'Found {len(sample_ids_not_in_mt)}/{len(sample_ids)} samples '
            f'in the subset set that do not matching IDs in the variant callset.\n'
            f'IDs that aren\'t in the callset: {sample_ids_not_in_mt}\n'
            f'All callset sample IDs: {mt_sample_ids}',
        )

    logging.info(
        f'Found {len(mt_sample_ids)} samples in mt, '
        f'subsetting to {len(sample_ids)} samples.'
    )

    n_rows_before = mt.count_rows()

    mt = mt.filter_cols(hl.literal(sample_ids).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logging.info(
        f'Finished subsetting to {len(sample_ids)} samples. '
        f'Kept {mt.count_cols()}/{len(mt_sample_ids)} samples, '
        f'{mt.count_rows()}/{n_rows_before} rows'
    )
    mt.write(str(out_mt_path), overwrite=True)
    logging.info(f'Written {out_mt_path}')


def annotate_dataset_mt(mt_path, out_mt_path, checkpoint_prefix):
    """
    Add dataset-level annotations.
    """
    mt = hl.read_matrix_table(str(mt_path))

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
        'dp': hl.if_else(
            is_called, hl.int(hl.min(mt.DP, 32000)), hl.missing(hl.tfloat)
        ),
        'sample_id': mt.s,
    }
    logging.info('Annotating genotypes')
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(hl.struct(**genotype_fields)),
    )

    mt = mt.checkpoint(f'{checkpoint_prefix}/dataset-genotypes.mt')
    logging.info(f'Written {checkpoint_prefix}/dataset-genotypes.mt')

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
        samples_num_alt=hl.struct(
            **{
                ('%i' % i): _genotype_filter_samples(_filter_num_alt(i))
                for i in range(1, 3, 1)
            }
        ),
        samples_gq=hl.struct(
            **{
                ('%i_to_%i' % (i, i + 5)): _genotype_filter_samples(
                    _filter_samples_gq(i)
                )
                for i in range(0, 95, 5)
            }
        ),
        samples_ab=hl.struct(
            **{
                '%i_to_%i' % (i, i + 5): _genotype_filter_samples(_filter_samples_ab(i))
                for i in range(0, 45, 5)
            }
        ),
    )
    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logging.info(f'Written {out_mt_path}')
