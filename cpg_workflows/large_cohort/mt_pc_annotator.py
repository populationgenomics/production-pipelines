#!/usr/bin/env python3

"""
analysis-runner --access-level "test" --dataset "bioheart" --description "QC PCA annotator" --output-dir "str/qc/iterative_pca_input" mt_pc_annotator.py \
--mt-path=gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt

"""

import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch,output_path

from ast import literal_eval
from cpg_utils import to_path


config = get_config()

@click.option('--mt-path', help='GCS Path to the input MT')
@click.command()
def main(mt_path):
    init_batch()

    mt = hl.read_matrix_table(mt_path)

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
    mt = mt.filter_cols(
        (mt.geno_pc1 >= -0.05) & (mt.geno_pc6 <= 0.05) & (mt.geno_pc6 >= -0.05)
    )

    with to_path(
        'gs://cpg-bioheart-test/str/associatr/input_files/remove-samples.txt'
    ).open() as f:
        array_string = f.read().strip()
        remove_samples = literal_eval(array_string)

    # remove related individuals
    mt = mt.filter_cols(hl.literal(remove_samples).contains(mt.s), keep=False)

    # drop chrX
    mt = mt.filter_rows((hl.str(mt.locus.contig).startswith('chrX')), keep=False)

    mt = mt.annotate_rows(motif_length = hl.len(mt.info.RU))
    mt = mt.annotate_rows(locus_length = mt.info.REF * mt.motif_length)

    #restrict to 2-6 bp motifs
    mt = mt.filter_rows((mt.motif_length >= 2) & (mt.motif_length <= 6))

    #remove segdup regions
    segdups = hl.import_bed('gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/segDupRegions/segDupRegions_hg38_sorted.bed.gz', force_bgz = True)
    mt = mt.annotate_rows(segdup_region = hl.is_defined(segdups[mt.locus]))
    mt = mt.filter_rows(mt.segdup_region == False)

    # tighten hwep
    annotations = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt')
    annotation_table = annotations.rows()
    mt = mt.annotate_rows(binom_hwep = annotation_table[mt.info.VARID].binom_hwep)
    mt = mt.filter_rows(mt.binom_hwep >= 0.05)

    # drop rows with locus length >100 bp
    mt = mt.filter_rows(mt.locus_length <= 100)


    # calculate the summed repeat length
    mt = mt.annotate_entries(sum_length=mt.allele_1_rep_length + mt.allele_2_rep_length)

    mt = mt.annotate_rows(mean_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.mean(mt.sum_length)),
                      stdev_sum_length = hl.agg.filter(hl.is_defined(mt.sum_length), hl.agg.stats(mt.sum_length)[1]))


    # drop rows with missing sum_length
    #missing_condition = hl.is_missing(mt.sum_length)
    #mt =mt.annotate_rows(missing_count=hl.agg.count_where(missing_condition))
    #mt = mt.filter_rows(mt.missing_count == 0)

    # write out mt to GCS path
    mt.write(output_path('str_pc_option7_annotated.mt'), overwrite=True)

    # print mt schema
    mt.describe()

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter


