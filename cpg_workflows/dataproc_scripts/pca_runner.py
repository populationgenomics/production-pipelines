import csv

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path


def pca_runner(vds_path, sample_id_file_path):
    vds = hl.vds.read_vds(str(vds_path))
    mt = vds.variant_data

    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
    if sample_id_file_path:
        with to_path(str(sample_id_file_path)).open('r') as f:
            sample_ids = list(csv.reader(f))[0]
        mt = mt.filter_cols(hl.literal(sample_ids).contains(mt.s))

    # subset to sites table
    sites_table = "gs://cpg-common-main/references/ancestry/pruned_variants.ht"
    qc_variants_ht = hl.read_table(sites_table)
    qc_variants_ht = qc_variants_ht.key_by('locus')

    mt = mt.filter_rows(hl.is_defined(qc_variants_ht[mt.locus]))

    mt = mt.checkpoint(output_path('checkpoint.mt', 'tmp'))
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=True)

    scores_output_path = output_path('scores.tsv', 'analysis')
    scores.export(scores_output_path)

    loadings_output_path = output_path('loadings.tsv', 'analysis')
    loadings.export(str(loadings_output_path))

    # Convert the list to a regular Python list
    eigenvalues_list = hl.eval(eigenvalues)
    # write the eigenvalues to a file
    with to_path(output_path('eigenvalues.txt', 'analysis')).open('w') as f:
        for item in eigenvalues_list:
            f.write(f'{item}\n')
