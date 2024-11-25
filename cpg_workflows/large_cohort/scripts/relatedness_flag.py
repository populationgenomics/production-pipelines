import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_config
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop


def cli_main():
    """
    A command-line entrypoint for the pcrelate process
    """

    parser = ArgumentParser()
    parser.add_argument('--relatedness', help='The Hail Table of sample relatedness')
    parser.add_argument('--qc', help='The Hail Table of Sample QC')
    parser.add_argument('--ht-out', help='The output path to write to')
    parser.add_argument('--gcta-out', help='The output path to write to')
    parser.add_argument('--checkpoint', help='Optional path to a Checkpointing directory', default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    flag_related(
        pcrelate_ht_path=args.relatedness,
        qc_ht_path=args.qc,
        ht_out=args.ht_out,
        gcta_out=args.gcta_out,
        checkpoint=args.checkpoint,
    )


def write_sample_ids_to_file(sample_ids: set, output_path: str):
    """
    GCTA PCA requires a list of id's to remove. Unable to read in a Hail Table in Batch so reading and writing here.
    Write a list of sample IDs to a file, one per line.
    """
    with to_path(output_path).open('w') as f:
        for sample_id in sample_ids:
            f.write(f'{sample_id}\n')


def flag_related(pcrelate_ht_path: str, qc_ht_path: str, ht_out: str, gcta_out: str, checkpoint: str | None = None):
    """
    Rank samples and flag samples to drop so there is only one sample per family
    left, with the highest rank in the family. The ranking is based on either
    'var_data_chr20_mean_dp' or 'autosomal_mean_dp' depending on which is present
    in the `sample_qc_ht` table.

    Args:
        pcrelate_ht_path (str): table with relatedness information.
        qc_ht_path (str): table with `var_data_chr20_mean_dp` or `autosomal_mean_dp` and `filters` fields.
        ht-out (str): path to write the output table of samples to drop.
        gcta-out (str): path to write the output table of samples to drop for gcta.
        checkpoint (str): optional, path to a checkpoint directory
    """
    logging.info('Flagging related samples to drop')

    ncpu = config_retrieve(['RelatednessFlag', 'cores'], 8)
    hl.context.init_spark(master=f'local[{ncpu}]', quiet=True)
    hl.default_reference('GRCh38')

    # read the sample HT
    sample_qc_ht = hl.read_table(qc_ht_path)

    # compute sample rankings
    rank_ht = _compute_sample_rankings(sample_qc_ht)

    # optionally, smash out a checkpoint
    if checkpoint:
        rank_ht = rank_ht.checkpoint(join(checkpoint, 'rank.ht'))

    # Remove only related individuals from the PCA, or related individuals AND
    # samples that did not meet sample QC thresholds
    if get_config()['large_cohort']['remove_failed_qc_pca']:
        try:
            filtered_samples = hl.literal(rank_ht.aggregate(hl.agg.filter(rank_ht.filtered, hl.agg.collect(rank_ht.s))))
        except hl.ExpressionException:
            # Hail doesn't handle it with `aggregate` when none of
            # the samples is 'filtered'
            filtered_samples = hl.empty_array(hl.tstr)
    else:
        filtered_samples = None

    # read relatedness data from a path
    relatedness_ht = hl.read_table(pcrelate_ht_path)

    to_drop_ht = compute_related_samples_to_drop(
        relatedness_ht,
        rank_ht,
        kin_threshold=get_config()['large_cohort']['max_kin'],
        filtered_samples=filtered_samples,
    )

    write_sample_ids_to_file(to_drop_ht.s.collect(), gcta_out)
    to_drop_ht.write(ht_out)


def _compute_sample_rankings(ht: hl.Table) -> hl.Table:
    """
    Orders samples by hard filters and coverage and adds rank, which is the lower,
    the better.

    @param ht: table with a `var_data_chr20_mean_dp` and `filters` fields.
    @return: table ordered by rank, with the following row fields:
        `rank`, `filtered`
    """
    ht = ht.drop(*list(ht.globals.dtype.keys()))
    ht = ht.select('var_data_chr20_mean_dp', filtered=hl.len(ht.filters) > 0)
    ht = ht.order_by(ht.filtered, hl.desc(ht.var_data_chr20_mean_dp))
    ht = ht.add_index(name='rank')
    return ht.key_by('s').select('filtered', 'rank')


if __name__ == '__main__':
    cli_main()
