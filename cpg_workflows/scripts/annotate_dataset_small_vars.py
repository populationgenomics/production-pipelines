"""
Read in a MT, and re-jig the annotations ready for Seqr Export
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.utils import get_logger


def annotate_dataset_mt(mt_path: str, out_mt_path: str):
    """
    Add dataset-level annotations.
    """

    init_batch()

    # this overrides the jar spec for the current session
    # and requires `init_batch()` to be called before any other hail methods
    if jar_spec := config_retrieve(['workflow', 'jar_spec_revisions', 'annotate_dataset'], False):
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
    get_logger().info('Annotating genotypes')
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
    get_logger().info(f'Written {out_mt_path}')


def cli_main():
    """
    CLI entrypoint
    """
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input MatrixTable to subset')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    args = parser.parse_args()
    annotate_dataset_mt(
        mt_path=args.input,
        out_mt_path=args.output,
    )


if __name__ == '__main__':
    cli_main()
