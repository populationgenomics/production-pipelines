import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils.config import config_retrieve


def cli_main():
    """
    A command-line entrypoint for the pcrelate process
    """

    parser = ArgumentParser()
    parser.add_argument('--dense_mt', help='The densified MatrixTable of input variants')
    parser.add_argument('--out', help='The output path to write to')
    parser.add_argument('--checkpoint', help='Optional path to a Checkpointing directory', default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    pcrelate(dense_mt_path=args.dense_mt, relatedness_out=args.out, checkpoint_dir=args.checkpoint)


def pcrelate(dense_mt_path: str, relatedness_out: str, checkpoint_dir: str | None = None):
    """
    Writes table with the following structure:
    Row fields:
        'i': str
        'j': str
        'kin': float64
        'ibd0': float64
        'ibd1': float64
        'ibd2': float64
    Key: ['i', 'j']

    Args:
        dense_mt_path ():
        relatedness_out ():
        checkpoint_dir (): a path to write a checkpoint to
    """

    # start up a Hail runtime
    ncpu = config_retrieve(['RelatednessPCRelate', 'cores'], 8)
    hl.context.init_spark(master=f'local[{ncpu}]', quiet=True)
    hl.default_reference('GRCh38')

    mt = hl.read_matrix_table(dense_mt_path).select_entries('GT')

    logging.info('Running relatedness check')

    sample_num = mt.cols().count()
    _, scores_ht, _ = hl.hwe_normalized_pca(mt.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False)

    # checkpoint if we were iven a location to write to
    if checkpoint_dir:
        checkpoint_path = join(checkpoint_dir, 'scores.ht')
        scores_ht.checkpoint(checkpoint_path, overwrite=True)

    relatedness_ht = hl.pc_relate(
        mt.GT,
        min_individual_maf=0.01,
        scores_expr=scores_ht[mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
    )

    # Converting keys for type struct{str} to str to align with the rank_ht `s` key
    relatedness_ht.key_by(i=relatedness_ht.i.s, j=relatedness_ht.j.s).write(relatedness_out)


if __name__ == '__main__':
    cli_main()
