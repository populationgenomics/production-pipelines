
import hail as hl

from cpg_utils import to_path


def pca_runner(file_path):

    mt = hl.read_matrix_table(str(file_path))

    #mt = mt.annotate_entries(
    #    allele_1_rep_length=hl.int(mt.REPCN.split('/')[0]),
    #    allele_2_rep_length=hl.if_else(
    #        hl.len(mt.REPCN.split('/')) == 2,
    #       hl.int(mt.REPCN.split('/')[1]),
    #        hl.missing('int32'),
    #    ),
    #)

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

    scores_output_path = 'gs://cpg-bioheart-test/str/qc/polymorphic_run_n2045/str_pca/scores.tsv.bgz'
    scores.export(str(scores_output_path))

    loadings_output_path = 'gs://cpg-bioheart-test/str/qc/polymorphic_run_n2045/str_pca/loadings.tsv.bgz'
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path('gs://cpg-bioheart-test/str/qc/polymorphic_run_n2045/str_pca/eigenvalues.txt').open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')