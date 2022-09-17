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
from hail_scripts.computed_fields import vep, variant_id

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='DEBUG', fmt=fmt)


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


def load_vqsr(site_only_vqsr_vcf_path: str):
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
        reference_genome=genome_build(),
    ).rows()

    # Annotating FILTERS from AS_FilterStatus before dropping that latter
    ht = ht.annotate(
        filters=ht.filters.union(hl.set(ht.info.AS_FilterStatus)).filter(
            lambda val: val != 'PASS'
        ),
    )

    # Dropping all INFO/AS* annotations as well as InbreedingCoeff, as they are
    # causing problems splitting multiallelics after parsing by Hail.
    as_fields = [f for f in ht.info if f.startswith('AS_')]
    fields_to_drop = as_fields + ['InbreedingCoeff']
    ht = ht.annotate(info=ht.info.drop(*fields_to_drop))

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

    ref_ht = hl.read_table(str(reference_path('seqr/combined_reference')))
    clinvar_ht = hl.read_table(str(reference_path('seqr/clinvar')))

    logging.info('Annotating with seqr-loader fields: round 1')
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
