"""
Hail query functions for seqr loader.
"""

import logging
import hail as hl
from lib.model.seqr_mt_schema import SeqrGenotypesSchema

logger = logging.getLogger(__file__)


class SeqrGenotypesESSchema(SeqrGenotypesSchema):
    """
    Modified version of SeqrVariantsAndGenotypesSchema to just base
    on SeqrGenotypesSchema, as we apply SeqrVariantSchema separately
    on the cohort level
    """

    @staticmethod
    def elasticsearch_row(ds):
        """
        Prepares the mt to export using ElasticsearchClient V02.
        - Flattens nested structs
        - drops locus and alleles key
        """
        # Converts a mt to the row equivalent.
        if isinstance(ds, hl.MatrixTable):
            ds = ds.rows()
        ds.describe()
        # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
        table = ds.drop('vep').flatten()
        # When flattening, the table is unkeyed, which causes problems because our locus and alleles should not
        # be normal fields. We can also re-key, but I believe this is computational?
        table = table.drop(table.locus, table.alleles)
        return table


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

    logger.info(
        f'Found {len(mt_sample_ids)} samples in mt, '
        f'subsetting to {len(sample_ids)} samples.'
    )

    n_rows_before = mt.count_rows()

    mt = mt.filter_cols(hl.literal(sample_ids).contains(mt.s))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info(
        f'Finished subsetting to {len(sample_ids)} samples. '
        f'Kept {mt.count_cols()}/{len(mt_sample_ids)} samples, '
        f'{mt.count_rows()}/{n_rows_before} rows'
    )
    mt.write(str(out_mt_path), overwrite=True)
    logger.info(f'Written {out_mt_path}')


def annotate_dataset_mt_old(mt_path, out_mt_path, checkpoints_bucket, overwrite=False):
    """
    Add dataset-level annotations.
    """
    mt = hl.read_matrix_table(str(mt_path))

    annotation_schema = SeqrGenotypesESSchema(mt)
    mt = annotation_schema.annotate_all(overwrite=True).mt

    mt = mt.checkpoint(
        f'{checkpoints_bucket}/dataset-genotypes.mt', _read_if_exists=not overwrite
    )
    logger.info(f'Written {checkpoints_bucket}/dataset-genotypes.mt')

    def _genotype_filter_samples(fn):
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    mt.describe()
    mt.write(out_mt_path, overwrite=True)
    logger.info(f'Written {out_mt_path}')
