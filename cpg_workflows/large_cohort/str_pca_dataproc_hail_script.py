
import hail as hl

from cpg_utils import to_path


def pca_runner(file_path):

    mt = hl.read_matrix_table(str(file_path))

    mt = mt.annotate_entries(
        allele_1_rep_length=hl.int(mt.REPCN.split('/')[0]),
        allele_2_rep_length=hl.if_else(
            hl.len(mt.REPCN.split('/')) == 2,
           hl.int(mt.REPCN.split('/')[1]),
            hl.missing('int32'),
        ),
    )

    # drop chrX
    mt = mt.filter_rows((hl.str(mt.locus.contig).startswith('chrX')), keep=False)

    mt = mt.annotate_rows(motif_length = hl.len(mt.info.RU))
    mt = mt.annotate_rows(locus_length = mt.info.REF * mt.motif_length)

    #restrict to 2-6 bp motifs
    #mt = mt.filter_rows((mt.motif_length >= 2) & (mt.motif_length <= 6))

    #remove segdup regions
    #segdups = hl.import_bed('gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/segDupRegions/segDupRegions_hg38_sorted.bed.gz', force_bgz = True)
    #mt = mt.annotate_rows(segdup_region = hl.is_defined(segdups[mt.locus]))
    #mt = mt.filter_rows(mt.segdup_region == False)

    # tighten hwep
    #annotations = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt')
    #annotation_table = annotations.rows()
    #mt = mt.annotate_rows(binom_hwep = annotation_table[mt.info.VARID].binom_hwep)
    #mt = mt.filter_rows(mt.binom_hwep >= 0.05)

    # drop rows with locus length >150 bp
    mt = mt.filter_rows(mt.locus_length <= 150)


    # calculate the summed repeat length
    mt = mt.annotate_entries(sum_length=mt.allele_1_rep_length + mt.allele_2_rep_length)

    mt = mt.annotate_rows(mean_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.mean(mt.sum_length)),
                      stdev_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.stats(mt.sum_length)[1]))


    # drop rows with missing sum_length
    missing_condition = hl.is_missing(mt.sum_length)
    mt =mt.annotate_rows(missing_count=hl.agg.count_where(missing_condition))
    mt = mt.filter_rows(mt.missing_count == 0)

    #  replace missing sum_length with the mean for that locus (PCA doesn't accept missing values) and normalise
    mt = mt.annotate_entries(
        sum_length_normalised=hl.or_else((mt.sum_length - mt.mean_sum_length) / mt.stdev_sum_length, 0.0))


    # run PCA
    eigenvalues, scores, loadings = hl.pca(mt.sum_length_normalised, k=10, compute_loadings=True)

    scores_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_5/scores.tsv.bgz'
    scores.export(str(scores_output_path))

    loadings_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_5/loadings.tsv.bgz'
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path('gs://cpg-bioheart-test/str/qc/iterative_pca/option_5/eigenvalues.txt').open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')