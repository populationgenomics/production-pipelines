#!/usr/bin/env python3
# pylint: skip-file

import hailtop.batch as hb
from os.path import join
import os


CONTAINER = f'australia-southeast1-docker.pkg.dev/cpg-common/images/alignment:v4'
REF_BUCKET = 'gs://cpg-reference/hg38/v1'
TARGET_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
TARGET_FASTA = join(TARGET_BUCKET, 'Homo_sapiens_assembly38.fasta')


def _index_bwa_job(b: hb.Batch, reference: hb.ResourceGroup):
    # exts = ['sa', 'amb', 'bwt', 'ann', 'pac']
    exts = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac']
    j = b.new_job('Index BWA')
    j.image(CONTAINER)
    total_cpu = 32
    j.cpu(total_cpu)
    j.storage('40G')
    j.declare_resource_group(bwa_index={e: '{root}.' + e for e in exts})
    j.command(
        f"""
set -o pipefail
set -ex

bwa-mem2 index {reference.base} -p {j.bwa_index}

df -h; pwd; ls | grep -v proc | xargs du -sh
    """
    )
    b.write_output(j.bwa_index, TARGET_FASTA)
    return j


def _test_bwa(
    b: hb.Batch,
    reference: hb.ResourceGroup,
):
    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test BWA')
    j.image(CONTAINER)
    total_cpu = 16
    bamsormadup_cpu = 3
    bwa_cpu = 1
    j.memory('highmem')
    j.cpu(total_cpu)
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.crai',
        }
    )
    sn = 'TEST'
    rg_line = f'@RG\\tID:{sn}\\tSM:~{sn}'
    # bwa index {reference.base}
    use_bazam = True

    sorted_bam = f'$(dirname {j.output_cram.cram})/sorted.bam'
    j.command(
        f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {fq1} - | \\
samtools sort -T $(dirname {j.output_cram.cram})/samtools-sort-tmp -Obam -o {sorted_bam}

picard MarkDuplicates I={sorted_bam} O=/dev/stdout M={j.duplicate_metrics} \\
    ASSUME_SORT_ORDER=coordinate | \\
samtools view -@30 -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
    """
    )
    return j


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ.get('HAIL_BUCKET', 'cpg-seqr-test-tmp')
print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='test')
reference = b.read_input_group(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
)
j1 = _index_bwa_job(b, reference)
j2 = _test_bwa(b, reference)
j2.depends_on(j1)
b.run(open=True)
