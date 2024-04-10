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

    annotations = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt')
    annotation_table = annotations.rows()
    mt = mt.annotate_rows(variant_lc = annotation_table[mt.info.VARID].variant_lc)
    mt = mt.annotate_rows(binom_hwep = annotation_table[mt.info.VARID].binom_hwep)
    mt = mt.annotate_rows(motif_length = annotation_table[mt.info.VARID].motif_length)
    mt = mt.annotate_rows(mode_allele = annotation_table[mt.info.VARID].aggregated_info.mode_allele)
    mt = mt.annotate_rows(prop_alleles_is_not_mode = annotation_table[mt.info.VARID].prop_alleles_is_not_mode)
    mt = mt.annotate_rows(num_alleles = annotation_table[mt.info.VARID].num_alleles)
    mt = mt.annotate_rows(obs_het = annotation_table[mt.info.VARID].obs_het)

    mt.rows().export('gs://cpg-bioheart-test/str/qc/iterative_pca_input/filtered_mt_rows.tsv')



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter


