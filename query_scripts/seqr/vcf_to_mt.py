#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Converts input multi-sample VCFs into a matrix table, and annotates it.
"""
from enum import Enum
from os.path import join
import click
import logging
import hail as hl
from gnomad.utils.sparse_mt import split_info_annotation
from lib.model.seqr_mt_schema import SeqrVariantSchema
from lib.model.base_mt_schema import row_annotation, RowAnnotationOmit

logger = logging.getLogger(__file__)

GENOME_VERSION = 'GRCh38'
REF_BUCKET = 'gs://cpg-reference'
BROAD_REF_BUCKET = f'{REF_BUCKET}/hg38/v1'
SEQR_REF_BUCKET = 'gs://cpg-seqr-reference-data'
REF_HT = f'{SEQR_REF_BUCKET}/GRCh38/all_reference_data/v2/combined_reference_data_grch38-2.0.3.ht'
CLINVAR_HT = f'{SEQR_REF_BUCKET}/GRCh38/clinvar/clinvar.GRCh38.2020-06-15.ht'
NAGIM_FREQS_HT = f'{REF_BUCKET}/seqr/v0-1-test/nagim-frequencies.ht'


SNP_SCORE_CUTOFF = 0
INDEL_SCORE_CUTOFF = 0


@click.command()
@click.option(
    '--vcf-path',
    'vcf_path',
    required=True,
)
@click.option(
    '--site-only-vqsr-vcf-path',
    'site_only_vqsr_vcf_path',
    required=True,
)
@click.option(
    '--dest-mt-path',
    'dest_path',
    required=True,
)
@click.option(
    '--bucket',
    'work_bucket',
    required=True,
)
@click.option('--disable-validation', 'disable_validation', is_flag=True)
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate matrix tables',
)
@click.option(
    '--overwrite',
    'overwrite',
    is_flag=True,
    help='Reuse intermediate files',
)
@click.option(
    '--run-vep',
    'run_vep',
    is_flag=True,
    help='Run VEP',
)
def main(
    vcf_path: str,
    site_only_vqsr_vcf_path: str,
    dest_path: str,
    work_bucket: str,
    disable_validation: bool,
    make_checkpoints: bool,
    overwrite: bool,
    run_vep: bool,
):  # pylint: disable=missing-function-docstring
    """
    Entry point
    """
    logger.info('Starting the seqr_load pipeline')

    hl.init(default_reference=GENOME_VERSION)

    mt = import_vcf(vcf_path, GENOME_VERSION)
    mt = annotate_old_and_split_multi_hts(mt)
    if not disable_validation:
        validate_mt(mt, sample_type='WGS')

    out_path = join(work_bucket, 'vqsr_and_37_coords.mt')
    if can_reuse(out_path, overwrite):
        mt = hl.read_matrix_table(out_path)
    else:
        vqsr_ht = load_vqsr(
            site_only_vqsr_vcf_path,
            output_ht_path=join(work_bucket, 'vqsr.ht'),
            overwrite=overwrite,
        )
        mt = annotate_vqsr(mt, vqsr_ht)
        mt = add_37_coordinates(mt)
        if make_checkpoints:
            mt.write(out_path, overwrite=True)
            mt = hl.read_matrix_table(out_path)

    if run_vep:
        out_path = join(work_bucket, 'vqsr_and_37_coords.vep.mt')
        if can_reuse(out_path, overwrite):
            mt = hl.read_matrix_table(out_path)
        else:
            mt = hl.vep(
                mt, 
                block_size=1000,
                # We are not starting the cluster with --vep, instead passing custom
                # startup script with --init gs://cpg-reference/vep/vep-GRCh38.sh,
                # so VEP_CONFIG_URI will not be set, thus need to provide config
                # as a function parameter here:
                # config='file:///vep_data/vep-gcloud.json'
            )
            if make_checkpoints:
                mt.write(out_path, overwrite=True)
                mt = hl.read_matrix_table(out_path)

    mt = compute_variant_annotated_vcf(
        mt, 
        ref_data=hl.read_table(REF_HT), 
        clinvar=hl.read_table(CLINVAR_HT), 
        nagim_freqs=hl.read_table(NAGIM_FREQS_HT),
    )
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=GENOME_VERSION.replace('GRCh', ''),
        sampleType='WGS',
        hail_version=hl.version(),
    )
    mt.describe()
    mt.write(dest_path, overwrite=True)


def annotate_vqsr(mt, vqsr_ht):
    """
    Assuming `vqsr_ht` is a site-only VCF from VQSR, and `mt` is a full matrix table,
    annotates `mt` rows from `vqsr_ht`
    """
    mt = mt.annotate_rows(**vqsr_ht[mt.row_key])
    
    # vqsr_ht has info annotation split by allele; plus new AS-VQSR annotations
    mt = mt.annotate_rows(info=vqsr_ht[mt.row_key].info)

    # populating filters which is outside of info
    mt = mt.annotate_rows(
        filters=mt.filters.union(vqsr_ht[mt.row_key].filters),
    )
    
    mt = mt.annotate_globals(**vqsr_ht.index_globals())
    return mt


def load_vqsr(
    site_only_vqsr_vcf_path: str,
    output_ht_path: str,
    overwrite: bool = False,
):
    """
    Loads the VQSR'ed site-only VCF into a site-only hail table. Populates ht.filters
    """
    if can_reuse(output_ht_path, overwrite):
        return hl.read_table(output_ht_path)

    logger.info(f'Importing VQSR annotations...')
    mt = hl.import_vcf(
        site_only_vqsr_vcf_path,
        force_bgz=True,
        reference_genome='GRCh38',
    ).repartition(1000)

    ht = mt.rows()

    # some numeric fields are loaded as strings, so converting them to ints and floats 
    ht = ht.annotate(
        info=ht.info.annotate(
            AS_VQSLOD=ht.info.AS_VQSLOD.map(hl.float),
            AS_QUALapprox=ht.info.AS_QUALapprox.split(r'\|')[1:].map(hl.int),
            AS_VarDP=ht.info.AS_VarDP.split(r'\|')[1:].map(hl.int),
            AS_SB_TABLE=ht.info.AS_SB_TABLE.split(r'\|').map(
                lambda x: hl.if_else(
                    x == '', hl.missing(hl.tarray(hl.tint32)), x.split(',').map(hl.int)
                )
            ),
        ),
    )
    unsplit_count = ht.count()

    ht = hl.split_multi_hts(ht)
    ht = ht.annotate(
        info=ht.info.annotate(**split_info_annotation(ht.info, ht.a_index)),
    )
    ht = ht.annotate(
        filters=ht.filters.union(hl.set([ht.info.AS_FilterStatus])),
    )
    ht.write(output_ht_path, overwrite=True)
    ht = hl.read_table(output_ht_path)
    logger.info(f'Wrote split HT to {output_ht_path}')
    split_count = ht.count()
    logger.info(
        f'Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations'
    )
    return ht


def import_vcf(source_paths, genome_version):
    """
    https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/lib/hail_tasks.py#L77-L89
    Import the VCFs from inputs. Set min partitions so that local pipeline execution
    takes advantage of all CPUs.
    :source_paths: list of paths to multi-sample VCFs
    :genome_version: GRCh37 or GRCh38
    @return a MatrixTable
    """
    recode = {}
    if genome_version == 'GRCh38':
        recode = {f'{i}': f'chr{i}' for i in (list(range(1, 23)) + ['X', 'Y'])}
    elif genome_version == 'GRCh37':
        recode = {f'chr{i}': f'{i}' for i in (list(range(1, 23)) + ['X', 'Y'])}

    return hl.import_vcf(
        source_paths,
        reference_genome=genome_version,
        skip_invalid_loci=True,
        contig_recoding=recode,
        force_bgz=True,
        min_partitions=500,
    )


def annotate_old_and_split_multi_hts(mt):
    """
    https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/seqr_loading.py#L89-L96

    Saves the old allele and locus because while split_multi does this, split_multi_hts drops this. Will see if
    we can add this to split_multi_hts and then this will be deprecated.
    @return: mt that has pre-annotations
    """
    # Named `locus_old` instead of `old_locus` because split_multi_hts drops `old_locus`.
    return hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )


class CPGSeqrVariantSchema(SeqrVariantSchema):
    """
    Modified schema to handle Nagim annotaion, and fix the AC annotation.
    """
    def __init__(self, *args, nagim_freqs: hl.Table, **kwargs):
        """
        Expects self._nagim_freqs to have the following fields:
        'freq': array<struct {
            AC: int32, 
            AF: float64, 
            AN: int32, 
        }> 
        'popmax': struct {
            AF: float64, 
        }
        """
        super().__init__(*args, **kwargs)
        ht = nagim_freqs.annotate(
            AF=nagim_freqs.freq.AF, 
            AC=nagim_freqs.freq.AC,
            AN=nagim_freqs.freq.AN,
            POPMAX_AF=nagim_freqs.popmax.AF,
        )
        ht = ht.select(ht.AF, ht.AC, ht.AN, ht.POPMAX_AF)
        self._nagim_data = ht

    @row_annotation(name='AC')
    def ac(self):  # pylint: disable=invalid-name,missing-function-docstring
        """
        AC in a split matrix table is a number, not an array.
        """
        return self.mt.info.AC

    @row_annotation
    def nagim(self):
        """
        Expects self._nagim_freqs to have: AF, AC, AN, POPMAX_AF
        """
        if self._nagim_data is None:
            raise RowAnnotationOmit

        return self._nagim_data[self.mt.row_key]


def compute_variant_annotated_vcf(
    mt: hl.MatrixTable,
    ref_data: hl.Table,
    clinvar: hl.Table,
    nagim_freqs: hl.Table,
) -> hl.MatrixTable:
    r"""
    Returns a matrix table with an annotated rows where each row annotation 
    is a previously called annotation (e.g. with the corresponding method or 
    all in case of `annotate_all`).

    Annotations are declared as methods on the schema_cls.

    There's a strong coupling between the @row_annotation decorator
    and the BaseMTSchema that's impractical to refactor, so we'll just leave it.

       class SomeClass(BaseMTSchema):
           @row_annotation()
           def a(self):
               return 'a_val'

    This loops through all the @row_annotation decorated methods
    on `schema_cls` and applies them all.
    
    See https://user-images.githubusercontent.com/22381693/113672350-f9b04000-96fa-11eb-91fe-e45d34c165c0.png
    for a rough overview of the structure and methods declared on:

                   BaseMTSchema
                     /       \
            SeqrSchema     SeqrGenotypesSchema
                    |         |
      SeqrVariantSchema       |
                    |         |
      CPGSeqrVariantSchema    |
                    \        /
            SeqrVariantsAndGenotypesSchema
            
    We apply only SeqrVariantSchema here (specifically, a modified 
    version CPGSeqrVariantSchema). SeqrGenotypesSchema is applied separately
    on the dataset level.
    """
    annotation_schema = CPGSeqrVariantSchema(
        mt, 
        ref_data=ref_data, 
        clinvar_data=clinvar, 
        nagim_freqs=nagim_freqs
    )
    mt = annotation_schema.annotate_all(overwrite=True).mt
    return mt


def elasticsearch_row(ds):
    """
    Copied from: https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/e41582d4842bc0d2e06b1d1e348eb071e00e01b3/luigi_pipeline/lib/model/seqr_mt_schema.py#L269-L290


    Prepares the mt to export using ElasticsearchClient V02.
    - Flattens nested structs
    - drops locus and alleles key

    TODO:
    - Call validate
    - when all annotations are here, whitelist fields to send instead of blacklisting.
    @return:
    """
    # Converts a mt to the row equivalent.
    if isinstance(ds, hl.MatrixTable):
        ds = ds.rows()
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    table = ds.drop('vep').flatten()
    # When flattening, the table is unkeyed, which causes problems because our locus and alleles should not
    # be normal fields. We can also re-key, but I believe this is computational?
    table = table.drop(table.locus, table.alleles)

    return table


def get_sample_type_stats(mt, threshold=0.3):
    """
    Calculate stats for sample type by checking against a list of common coding and non-coding variants.
    If the match for each respective type is over the threshold, we return a match.

    @param mt: Matrix Table to check
    @param threshold: if the matched percentage is over this threshold, we classify as match
    @return: a dict of coding/non-coding to dict with 'matched_count', 'total_count' and 'match' boolean.
    """
    stats = {}
    types_to_ht_path = {
        'noncoding': 'gs://seqr-reference-data/GRCh38/validate_ht/common_noncoding_variants.grch38.ht',
        'coding': 'gs://seqr-reference-data/GRCh38/validate_ht/common_coding_variants.grch38.ht',
    }
    for sample_type, ht_path in types_to_ht_path.items():
        ht = hl.read_table(ht_path)
        stats[sample_type] = ht_stats = {
            'matched_count': mt.semi_join_rows(ht).count_rows(),
            'total_count': ht.count(),
        }
        ht_stats['match'] = (
            ht_stats['matched_count'] / ht_stats['total_count']
        ) >= threshold
    return stats


class SeqrValidationError(Exception):
    """
    Thrown when the MatrixTable is failed
    """

    pass


def validate_mt(mt, sample_type):
    """
    Validate the mt by checking against a list of common coding and non-coding variants given its
    genome version. This validates genome_version, variants, and the reported sample type.

    @param mt: mt to validate
    @param sample_type: WGS or WES
    @return: True or Exception
    """
    sample_type_stats = get_sample_type_stats(mt)

    for name, stat in sample_type_stats.items():
        logger.info(
            'Table contains %i out of %i common %s variants.'
            % (stat['matched_count'], stat['total_count'], name)
        )

    has_coding = sample_type_stats['coding']['match']
    has_noncoding = sample_type_stats['noncoding']['match']

    if not has_coding and not has_noncoding:
        # No common variants detected.
        raise SeqrValidationError(
            'Genome version validation error: dataset specified but doesn\'t contain the expected number of common variants'
        )
    if has_noncoding and not has_coding:
        # Non coding only.
        raise SeqrValidationError(
            'Sample type validation error: Dataset contains noncoding variants but is missing common coding '
            'variants. Please verify that the dataset contains coding variants.'
        )
    if has_coding and not has_noncoding:
        # Only coding should be WES.
        if sample_type != 'WES':
            raise SeqrValidationError(
                f'Sample type validation error: dataset sample-type is specified as {sample_type} but appears to be '
                'WGS because it contains many common coding variants'
            )
    if has_noncoding and has_coding:
        # Both should be WGS.
        if sample_type != 'WGS':
            raise SeqrValidationError(
                f'Sample type validation error: dataset sample-type is specified as {sample_type} but appears to be '
                'WES because it contains many common non-coding variants'
            )

    return True


def add_37_coordinates(mt):
    """Annotates the GRCh38 MT with 37 coordinates using hail's built-in liftover
    @param mt: MatrixTable from VCF
    @return: MatrixTable annotated with GRCh37 coordinates
    """
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(join(REF_BUCKET, 'liftover/grch38_to_grch37.over.chain.gz'), rg37)
    mt = mt.annotate_rows(rg37_locus=mt.locus)
    return mt


class Cloud(Enum):
    """
    Cloud storage provider and correponding protocol prefix.
    """
    GS = 'gs'
    AZ = 'az'
    

def can_reuse(
    path,
    overwrite: bool,
    silent: bool = False,
) -> bool:
    return False


if __name__ == '__main__':
    main()  # pylint: disable=E1120
