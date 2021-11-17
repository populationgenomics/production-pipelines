#!/usr/bin/env python3

import hailtop.batch as hb
from os.path import join
import os

from cpg_pipes.jobs import wrap_command
from cpg_pipes.jobs.align import align, create_dragmap_index

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
DRAGMAP_IMAGE = f'{AR_REPO}/dragmap:1.2.1'
PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
DRAGMAP_INDEX_BUCKET = 'gs://cpg-reference/dragmap/v0'
DRAGMAP_INDEX_FILES = [
    'hash_table.cfg',
    'hash_table.cfg.bin',
    'hash_table.cmp',
    'hash_table_stats.txt',
    'ref_index.bin',
    'reference.bin',
    'repeat_mask.bin',
    'str_table.bin',
]


def _test_dragmap_job(b: hb.Batch):
    dragmap_index = b.read_input_group(
        **{
            k.replace('.', '_'): join(DRAGMAP_INDEX_BUCKET, k)
            for k in DRAGMAP_INDEX_FILES
        }
    )
    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test BWA')
    j.image(DRAGMAP_IMAGE)
    j.cpu(32)
    j.memory('standard')
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.crai',
        }
    )
    sn = 'TEST'
    j.command(wrap_command(f"""\
    dragen-os -r {dragmap_index} -1 {fq1} -2 {fq1} --RGID {sn} --RGSM {sn} |
    samtools sort -T $(dirname {j.sorted_bam})/samtools-sort-tmp -Obam -o {j.sorted_bam}
    """))
    return j


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ.get('HAIL_BUCKET', 'cpg-seqr-test-tmp')
print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
    token=os.environ.get('HAIL_TOKEN'),
)
b = hb.Batch(backend=backend, name='Create DRAGMAP index')
j1 = create_dragmap_index(b)
j2 = _test_dragmap_job(b)
j2.depends_on(j1)
b.run(wait=False)
