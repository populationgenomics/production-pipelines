"""
Use the gs://cpg-acute-care-test/mt/bbf37e78fe_5866-acute-care_full_copy.mt

Using the whole range of variants, this script will grab all the samples in the MT and do the following:

- generate a subset of 250 samples, export their data as a multisample VCF, and single-sample VCFs
- of that 250, generate a subset of 100, export as a multisample VCF
- repeat for 50, 25, and 10 samples

These multi/single sample data will be the basis of Nextflow benchmarking for the Talos annotation pre-process
"""
import logging
import os
import random

import hail as hl

from cpg_utils import hail_batch
from cpg_workflows import utils

logging.basicConfig(level=logging.INFO)

random.seed(42)

# kick up a batch
hail_batch.init_batch()

input_mt = 'gs://cpg-acute-care-test/mt/bbf37e78fe_5866-acute-care_full_copy.mt'
output_prefix = 'gs://cpg-acute-care-test/talos_benchmarking/ms_vcfs'
ss_vcf_prefix = 'gs://cpg-acute-care-test/talos_benchmarking/solo_vcfs'

mt = hl.read_matrix_table(input_mt)

samples = list(mt.s.collect())

# generate the biggest subset
_250_samples = random.sample(samples, 250)

# drop everything except variants
mt = mt.select_rows()

# select only the 250 samples
mt.filter_cols(hl.literal(_250_samples).contains(mt.s))

# strip out non-variant rows
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)
mt = mt.drop('variant_qc')

# write_this_little_mt, but seriously down-partitioned
mt = mt.repartition(50)
mt = mt.checkpoint('gs://cpg-acute-care-test/talos_benchmarking/250samples.mt', _read_if_exists=True)

for each_sam in _250_samples:
    out_path = os.path.join(ss_vcf_prefix, f'{each_sam}.vcf.bgz')

    if utils.exists(out_path):
        logging.info(f'{out_path} exists, skipping')
        continue

    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(mt.s == each_sam)
    hl.export_vcf(sam_mt, out_path, tabix=True)

# export the bigboi
out_path = os.path.join(output_prefix, '250.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_250_samples).contains(mt.s)), out_path, tabix=True)

# step down to 100 samples
_100_samples = random.sample(_250_samples, 100)
out_path = os.path.join(output_prefix, '100.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_100_samples).contains(mt.s)), out_path, tabix=True)

# step down to 50 samples
_50_samples = random.sample(_100_samples, 50)
out_path = os.path.join(output_prefix, '50.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_50_samples).contains(mt.s)), out_path, tabix=True)

# step down to 25 samples
_25_samples = random.sample(_50_samples, 25)
out_path = os.path.join(output_prefix, '25.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_25_samples).contains(mt.s)), out_path, tabix=True)

# step down to 10 samples
_10_samples = random.sample(_25_samples, 10)
out_path = os.path.join(output_prefix, '10.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_10_samples).contains(mt.s)), out_path, tabix=True)

# step down to 5 samples
_5_samples = random.sample(_10_samples, 5)
out_path = os.path.join(output_prefix, '5.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    hl.export_vcf(mt.filter_cols(hl.literal(_5_samples).contains(mt.s)), out_path, tabix=True)
