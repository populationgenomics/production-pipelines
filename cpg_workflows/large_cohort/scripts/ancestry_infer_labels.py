"""
This Stage is a thin wrapper around the gnomad methods PCA call
"""

import logging
import pickle
from argparse import ArgumentParser

import hail as hl

from cpg_utils.config import config_retrieve
from gnomad.sample_qc.ancestry import assign_population_pcs

MAX_MISLABELLED_TRAINING_SAMPLES: int = 50


def run_assign_population_pcs(pop_pca_scores_ht, min_prob, n_pcs: int):
    examples_num = pop_pca_scores_ht.aggregate(hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop)))
    logging.info(f'Running RF using {examples_num} training examples')
    pop_ht, pops_rf_model = assign_population_pcs(
        pop_pca_scores_ht,
        pc_cols=pop_pca_scores_ht.scores[:n_pcs],
        known_col='training_pop',
        min_prob=min_prob,
    )
    n_mislabeled_samples = pop_ht.aggregate(hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))
    return pop_ht, pops_rf_model, n_mislabeled_samples


def cli_main():
    """
    A command-line entrypoint for the ancestry add-background process
    """

    parser = ArgumentParser()
    parser.add_argument('--scores_ht', help='The PCA scores HT')
    parser.add_argument('--qc_in', help='The SampleQC HT to read')
    parser.add_argument('--qc_out', help='The updated SampleQC HT to write')
    parser.add_argument('--ht_out', help='population table output path')
    parser.add_argument('--pickle_out', help='population model output path')
    parser.add_argument('--txt_out', help='population inference output path')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(
        scores_ht_path=args.scores_ht,
        qc_in=args.qc_in,
        qc_out=args.qc_out,
        pop_ht_out=args.pop_ht_out,
        pickle_out=args.pickle_out,
        txt_out=args.txt_out,
    )


def main(scores_ht_path: str, qc_in: str, qc_out: str, pop_ht_out: str, pickle_out: str, txt_out: str):
    """

    Args:
        scores_ht_path ():
        qc_in (): where to read SampleQC HT from
        qc_out (): where to write updated SampleQC table
        pop_ht (): where to write the population table
        pickle_out (): where to write the model, pickled
        txt_out (): where to write population inference, text.gz
    """

    min_prob = config_retrieve(['large_cohort', 'min_pop_prob'])
    n_pcs = config_retrieve(['large_cohort', 'n_pcs'], 16)

    ncpu = config_retrieve(['Ancestry', 'cores'], 8)
    hl.context.init_spark(master=f'local[{ncpu}]', quiet=True)
    hl.default_reference('GRCh38')

    scores_ht = hl.read_table(scores_ht_path)
    sample_qc_ht = hl.read_table(qc_in)

    training_pop_ht = sample_qc_ht.filter(hl.is_defined(sample_qc_ht['training_pop']))

    # no idea if continuing is a good idea in this case - failing instead
    if training_pop_ht.count() < 2:
        raise ValueError('Need at least 2 samples with known `population` label to run PCA')

    logging.info(
        'Using calculated PCA scores as well as training samples with known '
        '`population` label to assign population labels to remaining samples',
    )
    scores_ht = scores_ht.annotate(training_pop=training_pop_ht[scores_ht.key].training_pop)

    pop_ht, pops_rf_model, n_mislabeled_samples = run_assign_population_pcs(scores_ht, min_prob, n_pcs)
    while n_mislabeled_samples > MAX_MISLABELLED_TRAINING_SAMPLES:
        logging.info(
            f'Found {n_mislabeled_samples} samples '
            f'labeled differently from their known pop. '
            f'Re-running without them.',
        )

        pop_ht = pop_ht[scores_ht.key]

        # not sure if persist is viable here
        pop_pca_scores_ht = scores_ht.annotate(
            training_pop=hl.or_missing((pop_ht.training_pop == pop_ht.pop), scores_ht.training_pop),
        ).persist()

        pop_ht, pops_rf_model, n_mislabeled_samples = run_assign_population_pcs(pop_pca_scores_ht, min_prob, n_pcs)

    # Writing a tab delimited file indicating inferred sample populations
    pc_cnt = min(hl.min(10, hl.len(pop_ht.pca_scores)).collect())
    pop_ht.transmute(**{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(pc_cnt)}).export(txt_out)

    # Writing the RF model used for inferring sample populations
    with open(pickle_out, 'wb', encoding='utf-8') as handle:
        pickle.dump(pops_rf_model, handle)

    pop_ht.annotate(is_training=hl.is_defined(training_pop_ht[pop_ht.key])).write(pop_ht_out)
    # read this version back in
    pop_ht = hl.read_table(pop_ht_out)

    # update contents and write out
    sample_qc_ht.annotate(**pop_ht[sample_qc_ht.key]).write(qc_out)


if __name__ == '__main__':
    cli_main()
