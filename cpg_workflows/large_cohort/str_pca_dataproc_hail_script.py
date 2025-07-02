from ast import literal_eval

import hail as hl

from cpg_utils import to_path
import pandas as pd


def pca_runner(file_path):
    mt = hl.read_matrix_table(str(file_path))

    #bioheart_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/bioheart_n975_sample_covariates.csv')['sample_id']
    ##tob_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/covariates/6_rna_pcs/CD4_TCM_covariates.csv')['sample_id']
    #samples = tob_ids.to_list()

    # filter the MT to only include samples in the sample list
    #mt = mt.filter_cols(hl.literal(samples).contains(mt.s))

    # remove monomorphic variants, set locus level call rate >=0.9, observed heterozygosity >=0.00995, locus level HWEP (binom definition) >=10^-6
    mt = mt.filter_rows(
        (mt.num_alleles > 1) & (mt.variant_qc.call_rate >= 0.9) & (mt.obs_het >= 0.00995) & (mt.binom_hwep >= 0.000001),
    )
    # set sample level call rate >=0.99
    mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.99)

    # set big expansions/deletions beyond [-30,20] relative to mode allele to NA
    condition_allele_1 = (mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)
    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_1,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_1,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.SO,
        ),
    )
    condition_allele_2 = (mt.allele_2_minus_mode > 20) | (mt.allele_2_minus_mode < -30)

    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_2,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_2,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.SO,
        ),
    )

    # calculate proportion of GTs that are defined per locus (after applying call-level filters, variant_qc.call_rate is not accurate anymore)
    mt = mt.annotate_rows(
        prop_GT_exists=hl.agg.count_where(hl.is_defined(mt.GT)) / (mt.variant_qc.n_called + mt.variant_qc.n_not_called),
    )
    # re-enforce locus level call rate >=0.9
    mt = mt.filter_rows(mt.prop_GT_exists >= 0.9)

    # create DS entry field (summed repeat length across both alleles); set to mode*2 if GT is missing
    mt = mt.annotate_entries(
        DS=hl.if_else(
            hl.is_defined(mt.GT),
            mt.allele_1_rep_length + mt.allele_2_rep_length,
            mt.aggregated_info.mode_allele * 2,
        ),
    )

    #mt = mt.annotate_rows(min_DS=hl.agg.min(mt.DS), max_DS=hl.agg.max(mt.DS))
    #mt = mt.annotate_entries(DS_rescaled=(mt.DS - mt.min_DS) / (mt.max_DS - mt.min_DS) * 2)
    # clean up mt (drop unnecessary fields)
    #mt = mt.drop('DS', 'min_DS', 'max_DS')

    # rename DS_rescaled as DS
    #mt = mt.annotate_entries(DS=mt.DS_rescaled)
    #mt = mt.drop('DS_rescaled')

    # table_geno_pcs = hl.import_table(
    #    'gs://cpg-bioheart-test/str/anndata/saige-qtl/input_files/covariates/sex_age_geno_pcs_tob_bioheart.csv',
    ##    delimiter=',',
    #    impute=True,
    # )

    # table_geno_pcs = table_geno_pcs.key_by('sample_id')
    # mt = mt.annotate_cols(sample_id = 'CPG'+hl.str(mt.s))
    # mt = mt.key_cols_by('sample_id')
    # mt = mt.annotate_cols(geno_pc1=hl.float(table_geno_pcs[mt.sample_id].geno_PC1))
    # mt = mt.annotate_cols(geno_pc6=hl.float(table_geno_pcs[mt.sample_id].geno_PC6))
    # remove ancestry outliers
    # mt = mt.filter_cols(
    #    (mt.geno_pc1 >=0.01) & (mt.geno_pc6 <= 0.04) & (mt.geno_pc6 >= -0.03)& (mt.geno_pc6 <= 0.01)
    # )

    with to_path('gs://cpg-bioheart-test/str/associatr/input_files/remove-samples.txt').open() as f:
        array_string = f.read().strip()
        remove_samples = literal_eval(array_string)

    # remove related individuals
    mt = mt.filter_cols(hl.literal(remove_samples).contains(mt.s), keep=False)

    # remove outlier
    # ids_to_filter = ['CPG309245', 'CPG315648','CPG312819','CPG316182','CPG311522','CPG315689','CPG315655','CPG310078']
    # mt = mt.filter_cols(hl.literal(ids_to_filter).contains(mt.s), keep = False)

    # remove putative variants driving the batch effect
    # table_variants = hl.import_table('gs://cpg-bioheart-test/str/filtered_variants_opt13.csv')
    # table_variants = table_variants.annotate(locus = hl.parse_locus(table_variants['locus']))
    # table_variants = table_variants.key_by('locus')

    # mt = mt.annotate_rows(not_batch_effect = hl.is_defined(table_variants[mt.locus]))
    # mt = mt.filter_rows(mt.not_batch_effect == True)

    # retain variants that pass min LC
    # table_gc = hl.import_table('gs://cpg-bioheart-test/str/qc/iterative_pca_input/REPID_passing_min_LC_7.5.csv')
    # table_gc = table_gc.key_by('REPID')
    # mt = mt.annotate_rows(passing_GC = hl.is_defined(table_gc[mt.info.REPID]))
    # mt = mt.filter_rows(mt.passing_GC == True)

    # drop chrX
    mt = mt.filter_rows((hl.str(mt.locus.contig).startswith('chrX')), keep=False)

    # mt = mt.annotate_rows(motif_length = hl.len(mt.info.RU))
    # mt = mt.annotate_rows(locus_length = mt.info.REF * mt.motif_length)

    # restrict to 2-6 bp motifs
    # mt = mt.filter_rows((mt.motif_length >= 2) & (mt.motif_length <= 6))

    # remove segdup regions
    # segdups = hl.import_bed('gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/segDupRegions/segDupRegions_hg38_sorted.bed.gz', force_bgz = True)
    # mt = mt.annotate_rows(segdup_region = hl.is_defined(segdups[mt.locus]))
    # mt = mt.filter_rows(mt.segdup_region == False)

    # only keep STRs in the EnsembleTR catalog intervals
    # illumina = hl.import_bed('gs://cpg-bioheart-test/str/batch_debug/EnsembleTR.bed')
    # mt = mt.annotate_rows(illumina_region = hl.is_defined(illumina[mt.locus]))
    # mt = mt.filter_rows(mt.illumina_region == True)

    # tighten hwep
    # annotations = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt')
    # annotation_table = annotations.rows()
    # mt = mt.annotate_rows(binom_hwep = annotation_table[mt.info.VARID].binom_hwep)
    # mt = mt.filter_rows(mt.binom_hwep >= 0.05)

    # drop rows with locus length >100 bp
    # mt = mt.filter_rows(mt.locus_length <= 100)

    # calculate the summed repeat length
    # mt = mt.annotate_entries(sum_length=mt.allele_1_rep_length + mt.allele_2_rep_length)

    # mt = mt.annotate_rows(mean_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.mean(mt.sum_length)),
    # stdev_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.stats(mt.sum_length)[1]))
    mt = mt.annotate_rows(
        mean_DS=hl.agg.filter(hl.is_defined(mt.DS), hl.agg.mean(mt.DS)),
        stdev_DS=hl.agg.filter(hl.is_defined(mt.DS), hl.agg.stats(mt.DS)[1]),
    )

    # remove non-variable loci
    mt = mt.filter_rows(mt.stdev_DS != 0)

    #  replace missing sum_length with the mean for that locus (PCA doesn't accept missing values) and normalise
    mt = mt.annotate_entries(DS_normalised=hl.or_else((mt.DS - mt.mean_DS) / mt.stdev_DS, 0.0))

    # run PCA
    eigenvalues, scores, loadings = hl.pca(mt.DS_normalised, k=10, compute_loadings=True)

    scores_output_path = 'gs://cpg-bioheart-test/str/pca/n950-tob-default-filters/scores.tsv.bgz'
    scores.export(str(scores_output_path))

    loadings_output_path = 'gs://cpg-bioheart-test/str/pca/n950-tob-default-filters/loadings.tsv.bgz'
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path('gs://cpg-bioheart-test/str/pca/n950-tob-default-filters/eigenvalues.txt').open(
        'w',
    ) as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')
