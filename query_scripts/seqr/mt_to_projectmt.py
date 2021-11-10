#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Subset matrix table to a list of samples and annotate with SeqrGenotypesSchema.
"""

import logging
import click
import hail as hl
from lib.model.seqr_mt_schema import SeqrGenotypesSchema

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--out-mt-path',
    'out_mt_path',
    required=True,
)
@click.option(
    '--subset-tsv',
    'subset_tsv_path',
    help='Path to a TSV file with one column of sample IDs: s.',
)
def main(
    mt_path: str,
    out_mt_path: str,
    subset_tsv_path: str,
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(mt_path)
    logger.info('Subsetting to the requested set of samples')
    if subset_tsv_path:
        mt = _subset_samples_and_variants(mt, subset_tsv_path)

    logger.info('Annotating genotypes')
    mt = _compute_genotype_annotated_vcf(mt)

    # Fixing NaN values in AS_MQ
    mt.annotate_rows(
        info=hl.struct(
            AS_MQ=hl.if_else(
                ~hl.is_nan(mt.info.AS_MQ),  # pylint: disable=invalid-unary-operand-type
                mt.info.AS_MQ,
                0,
            ),
        )
    )

    mt.describe()
    mt.write(out_mt_path, overwrite=True)


def _subset_samples_and_variants(mt, subset_tsv_path: str) -> hl.MatrixTable:
    """
    Subset the MatrixTable to the provided list of samples and to variants present
    in those samples
    :param mt: MatrixTable from VCF
    :param subset_tsv_path: path to a TSV file with a single column 's', with no header
    :return: MatrixTable subsetted to list of samples
    """
    count_before = mt.count_cols()

    subset_ht = hl.import_table(subset_tsv_path, no_header=True)
    subset_ht = subset_ht.transmute(s=subset_ht.f0).key_by('s')
    subset_count = subset_ht.count()
    anti_join_ht = subset_ht.anti_join(mt.cols())
    anti_join_ht_count = anti_join_ht.count()

    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        raise Exception(
            f'Only {subset_count - anti_join_ht_count} out of {count_before} '
            'subsetting-table IDs matched IDs in the variant callset.\n'
            f'IDs that aren\'t in the callset: {missing_samples}\n'
            f'All callset sample IDs:{mt.s.collect()}',
            missing_samples,
        )

    mt = mt.semi_join_cols(subset_ht)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info(
        f'Finished subsetting to {subset_count} samples. Kept {mt.count_cols()} '
        f'out of {count_before} samples in vds'
    )
    return mt


def _compute_genotype_annotated_vcf(
    mt, schema_cls=SeqrGenotypesSchema
) -> hl.MatrixTable:
    r"""
                    BaseMTSchema
                     /        \
            SeqrSchema         |
                    |          |
      SeqrVariantSchema     SeqrGenotypesSchema
                    |          |
      SeqrVariantASSchema     /
                    \        /
            SeqrVariantsAndGenotypesSchema

    SeqrVariantASSchema is applied on the cohort level separately.
    """
    annotation_schema = schema_cls(mt)
    mt = annotation_schema.annotate_all(overwrite=True).mt
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
