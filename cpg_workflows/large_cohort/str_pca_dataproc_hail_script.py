
import hail as hl
from ast import literal_eval


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
    table_geno_pcs = hl.import_table(
        'gs://cpg-bioheart-test/str/anndata/saige-qtl/input_files/covariates/sex_age_geno_pcs_tob_bioheart.csv',
        delimiter=',',
        impute=True,
    )

    table_geno_pcs = table_geno_pcs.key_by('sample_id')
    mt = mt.annotate_cols(sample_id = 'CPG'+hl.str(mt.s))
    mt = mt.key_cols_by('sample_id')
    mt = mt.annotate_cols(geno_pc1=hl.float(table_geno_pcs[mt.sample_id].geno_PC1))
    mt = mt.annotate_cols(geno_pc6=hl.float(table_geno_pcs[mt.sample_id].geno_PC6))
    # remove ancestry outliers
    #mt = mt.filter_cols(
    #    (mt.geno_pc1 >=0.01) & (mt.geno_pc6 <= 0.04) & (mt.geno_pc6 >= -0.03)& (mt.geno_pc6 <= 0.01)
    #)

    with to_path(
        'gs://cpg-bioheart-test/str/associatr/input_files/remove-samples.txt'
    ).open() as f:
        array_string = f.read().strip()
        remove_samples = literal_eval(array_string)

    # remove related individuals
    mt = mt.filter_cols(hl.literal(remove_samples).contains(mt.s), keep=False)

    # remove putative variants driving the batch effect
    table_variants = hl.import_table('gs://cpg-bioheart-test/str/filtered_variants_opt11.csv')
    table_variants = table_variants.annotate(locus = hl.parse_locus(table_variants['locus']))
    table_variants = table_variants.key_by('locus')

    mt = mt.annotate_rows(not_batch_effect = hl.is_defined(table_variants[mt.locus]))
    mt = mt.filter_rows(mt.not_batch_effect == True)


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

    # drop rows with locus length >100 bp
    #mt = mt.filter_rows(mt.locus_length <= 100)


    # calculate the summed repeat length
    mt = mt.annotate_entries(sum_length=mt.allele_1_rep_length + mt.allele_2_rep_length)

    mt = mt.annotate_rows(mean_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.mean(mt.sum_length)),
                      stdev_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.stats(mt.sum_length)[1]))

    #remove non-variable loci
    mt = mt.filter_rows(mt.stdev_sum_length != 0)

    #  replace missing sum_length with the mean for that locus (PCA doesn't accept missing values) and normalise
    mt = mt.annotate_entries(
        sum_length_normalised=hl.or_else((mt.sum_length - mt.mean_sum_length) / mt.stdev_sum_length, 0.0))

    # run PCA
    eigenvalues, scores, loadings = hl.pca(mt.sum_length_normalised, k=10, compute_loadings=True)

    scores_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_11/scores.tsv.bgz'
    scores.export(str(scores_output_path))

    loadings_output_path = 'gs://cpg-bioheart-test/str/qc/iterative_pca/option_11/loadings.tsv.bgz'
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path('gs://cpg-bioheart-test/str/qc/iterative_pca/option_11/eigenvalues.txt').open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')