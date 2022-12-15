#!/usr/bin/env python3

"""
Subset reference data to test intervals.
"""

import hail as hl
from cpg_utils import to_path
from cpg_utils.workflows.utils import exists
from cpg_utils.hail_batch import start_query_context


start_query_context(
    query_backend='batch',
    dataset='thousand-genomes',
    billing_project='thousand-genomes',
)

src_prefix = to_path('gs://cpg-common-main/references')
dst_prefix = to_path('gs://cpg-common-main/references/subset-toy-chr20-X-Y')

intervals_path = dst_prefix / 'intervals.bed'
if intervals_path.exists():
    intervals_path.unlink()
with (to_path(__file__) / 'intervals.bed').open() as f, intervals_path.open('w') as out:
    out.write(f.read())
intervals_ht = hl.import_bed(str(intervals_path))

for path in [
    'gnomad/v0/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht',
    'gnomad/v0/lcr_intervals/LCRFromHengHg38.ht',
    'gnomad/v0/seg_dup_intervals/GRCh38_segdups.ht',
    'gnomad/v0/clinvar/clinvar_20190923.ht',
    'gnomad/v0/hapmap/hapmap_3.3.hg38.ht',
    'gnomad/v0/kgp/1000G_omni2.5.hg38.ht',
    'gnomad/v0/kgp/1000G_phase1.snps.high_confidence.hg38.ht',
    'gnomad/v0/mills/Mills_and_1000G_gold_standard.indels.hg38.ht',
    'gnomad/v0/sample_qc/pre_ld_pruning_qc_variants.ht',
]:
    src_path = src_prefix / path
    dst_path = dst_prefix / path

    if exists(dst_path):
        continue
    if not src_path.suffix == '.ht':
        raise ValueError(f'Unrecognised extension: {src_path}')
    ht = hl.read_table(str(src_path))
    ht.describe()
    # noinspection PyProtectedMember
    if 'interval' in ht._fields:
        ht = ht.filter(hl.is_defined(intervals_ht[ht.interval]))
    else:
        # noinspection PyProtectedMember
        assert 'locus' in ht._fields
        ht = ht.filter(hl.is_defined(intervals_ht[ht.locus]))
    ht = ht.naive_coalesce(1)
    ht.write(str(dst_path), overwrite=True)
