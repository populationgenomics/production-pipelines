
import hail as hl
from ast import literal_eval


from cpg_utils import to_path


def pca_runner(file_path):

    mt = hl.read_matrix_table(str(file_path))

    # drop rows with missing sum_length
    missing_condition = hl.is_missing(mt.sum_length)
    mt =mt.annotate_rows(missing_count=hl.agg.count_where(missing_condition))
    mt = mt.filter_rows(mt.missing_count == 0)

    #  replace missing sum_length with the mean for that locus (PCA doesn't accept missing values) and normalise
    mt = mt.annotate_entries(
        sum_length_normalised=hl.or_else((mt.sum_length - mt.mean_sum_length) / mt.stdev_sum_length, 0.0))

    # run PCA
    eigenvalues, scores, loadings = hl._blanczos_pca(mt.sum_length_normalised, k=2, compute_loadings=True)

    scores_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_7/scores.tsv.bgz'
    scores.export(str(scores_output_path))

    loadings_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_7/loadings.tsv.bgz'
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path('gs://cpg-bioheart-test/str/qc/iterative_pca/option_7/eigenvalues.txt').open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')