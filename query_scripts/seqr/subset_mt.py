#!/usr/bin/env python3

"""
Hail Query script.

Subsets matrix table to a list of samples and annotate with SeqrGenotypesSchema.
"""

import logging
import click
import hail as hl
from lib.model.base_mt_schema import row_annotation
from lib.model.seqr_mt_schema import BaseMTSchema
from cpg_utils.hail import init_query_service

logger = logging.getLogger(__file__)


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
    """
    Entry point
    """
    init_query_service()

    mt = hl.read_matrix_table(mt_path)
    logger.info(f'Size of mt: {mt.count()}')
    
    logger.info('Subsetting to the requested set of samples')
    if subset_tsv_path:
        mt = _subset_samples_and_variants(mt, subset_tsv_path)
        logger.info(f'Size of mt after subsetting: {mt.count()}')

    logger.info('Annotating genotypes')
    mt = _compute_genotype_annotated_vcf(mt)

    # AS_MQ and InbreedingCoeff can be NaN, need to fix that to avoid ES loader 
    # failures like: 
    # `is.hail.relocated.org.elasticsearch.hadoop.rest.EsHadoopRemoteException: 
    # mapper_parsing_exception: failed to parse field [info_AS_MQ] of type [double]`
    # in: https://batch.hail.populationgenomics.org.au/batches/6621/jobs/6
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AS_MQ=hl.if_else(
                hl.is_nan(mt.info.AS_MQ), 
                0.0, 
                mt.info.AS_MQ
            ),
            InbreedingCoeff=hl.if_else(
                hl.is_nan(mt.info.InbreedingCoeff), 
                0.0, 
                mt.info.InbreedingCoeff
            ),
            # https://batch.hail.populationgenomics.org.au/batches/6973/jobs/12
            AS_InbreedingCoeff=hl.if_else(
                hl.is_nan(mt.info.AS_InbreedingCoeff), 
                0.0, 
                mt.info.AS_InbreedingCoeff
            ),
        )
    )

    mt.describe()
    mt.write(out_mt_path, overwrite=True)


def _subset_samples_and_variants(mt, subset_tsv_path: str) -> hl.MatrixTable:
    """
    Subset the MatrixTable to the provided list of samples and to variants present
    in those samples
    @param mt: MatrixTable from VCF
    @param subset_tsv_path: path to a TSV file with a single column 's', with no header
    @return: MatrixTable subsetted to list of samples
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


class SeqrGenotypesSchema(BaseMTSchema):
    @row_annotation()
    def genotypes(self):
        return hl.agg.collect(hl.struct(**self._genotype_fields()))

    @row_annotation(fn_require=genotypes)
    def samples_no_call(self):
        return self._genotype_filter_samples(lambda g: g.num_alt == -1)

    @row_annotation(fn_require=genotypes)
    def samples_num_alt(self, start=1, end=3, step=1):
        return hl.struct(
            **{
                ('%i' % i): 
                self._genotype_filter_samples(lambda g: g.num_alt == i)
                for i in range(start, end, step)
            }
        )

    # @row_annotation(fn_require=genotypes)
    # def samples_gq(self, start=0, end=95, step=5):
    #     # struct of x_to_y to a set of samples in range of x and y for gq.
    #     return hl.struct(
    #         **{
    #             ('%i_to_%i' % (i, i + step)): 
    #             self._genotype_filter_samples(
    #                 lambda g: ((g.gq >= i) & (g.gq < i + step))
    #             )
    #             for i in range(start, end, step)
    #         }
    #     )

    @row_annotation(fn_require=genotypes)
    def samples_ab(self, start=0, end=45, step=5):
        # struct of x_to_y to a set of samples in range of x and y for ab.
        return hl.struct(
            **{
                '%i_to_%i'
                % (i, i + step): self._genotype_filter_samples(
                    lambda g: (
                        (g.num_alt == 1)
                        & ((g.ab * 100) >= i)
                        & ((g.ab * 100) < i + step)
                    )
                )
                for i in range(start, end, step)
            }
        )

    def _genotype_filter_samples(self, filter):
        # Filter on the genotypes.
        return hl.set(self.mt.genotypes.filter(filter).map(lambda g: g.sample_id))

    def _genotype_fields(self):
        # Convert the mt genotype entries into num_alt, gq, ab, dp, and sample_id.
        is_called = hl.is_defined(self.mt.GT)
        return {
            'num_alt': hl.if_else(is_called, self.mt.GT.n_alt_alleles(), -1),
            'gq': hl.if_else(is_called, self.mt.GQ, hl.null(hl.tint)),
            'ab': hl.bind(
                lambda total: hl.if_else(
                    (is_called) & (total != 0) & (hl.len(self.mt.AD) > 1),
                    hl.float(self.mt.AD[1] / total),
                    hl.null(hl.tfloat),
                ),
                hl.sum(self.mt.AD),
            ),
            'dp': hl.if_else(
                is_called, hl.int(hl.min(self.mt.DP, 32000)), hl.null(hl.tfloat)
            ),
            'sample_id': self.mt.s,
        }


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
      CPGSeqrVariantSchema    /
                    \        /
            SeqrVariantsAndGenotypesSchema

    CPGSeqrVariantSchema is applied on the cohort level separately.
    """
    annotation_schema = schema_cls(mt)
    mt = annotation_schema.annotate_all(overwrite=True).mt
    return mt


if __name__ == '__main__':
    main()  # pylint: disable=E1120
