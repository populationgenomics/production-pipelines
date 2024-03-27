#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
Hail script to submit to a dataproc cluster

 analysis-runner --dataset "bioheart" \
    --description "str_pca" \
    --access-level "test" \
    --output-dir "str/qc/filtered_mt" \
    str_pca.py --file-path=gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt

"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch
from cpg_utils import to_path


@click.option(
    '--file-path',
    help='GCS file path to Hail STR matrix Table',
    type=str,
)
@click.command()
def pca_runner(file_path):

    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table(file_path)

    mt = mt.annotate_entries(
        allele_1_rep_length=hl.int(mt.REPCN.split('/')[0]),
        allele_2_rep_length=hl.if_else(
            hl.len(mt.REPCN.split('/')) == 2,
            hl.int(mt.REPCN.split('/')[1]),
            hl.missing('int32'),
        ),
    )

    # calculate the summed repeat length
    mt = mt.annotate_entries(sum_length=mt.allele_1_rep_length + mt.allele_2_rep_length)
    mt = mt.annotate_cols(mean_sum_length=hl.agg.mean(hl.float(mt.sum_length)))

    # replace missing sum_length with the mean for that locus (PCA doesn't accept missing values)
    mt = mt.annotate_entries(
        sum_length=hl.if_else(
            hl.is_missing(mt.sum_length), mt.mean_sum_length, mt.sum_length
        )
    )

    # mean-centre and normalise sum_length
    mt = mt.annotate_cols(sd_sum_length=hl.agg.stats(mt.sum_length).stdev)
    mt = mt.annotate_entries(
        sum_length=(mt.sum_length - mt.mean_sum_length) / mt.sd_sum_length
    )

    # run PCA
    eigenvalues, scores, loadings = hl.pca(mt.sum_length, k=10, compute_loadings=True)

    scores_output_path = output_path(f'str_pca/scores.tsv.bgz')
    scores.export(scores_output_path)

    loadings_output_path = output_path(f'str_pca/loadings.tsv.bgz')
    loadings.export(loadings_output_path)

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path(output_path(f'str_pca/eigenvalues.txt')).open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
