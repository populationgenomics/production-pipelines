#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config

import logging


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--sample-list',
    'sample_list_path',
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
    mt_path: str,
    sample_list_path: str,
    out_mt_path: str,
    checkpoint_prefix: str,
):
    hl.init(default_reference='GRCh38')

    with to_path(sample_list_path).open() as f:
        sample_ids = f.read().strip().split(',')

    subset_mt_path = to_path(checkpoint_prefix) / 'cohort-subset.mt'

    subset_mt_to_samples(
        mt_path=mt_path,
        sample_ids=sample_ids,
        out_mt_path=str(subset_mt_path),
    )

    annotate_dataset_mt(
        mt_path=str(subset_mt_path),
        out_mt_path=out_mt_path,
        checkpoint_prefix=checkpoint_prefix,
        overwrite=not get_config()['workflow'].get('check_intermediates'),
    )


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


def annotate_dataset_mt(mt_path, out_mt_path, checkpoint_prefix, overwrite=False):
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
    mt = mt.checkpoint(
        f'{checkpoint_prefix}/dataset-genotypes.mt', _read_if_exists=not overwrite
    )
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
