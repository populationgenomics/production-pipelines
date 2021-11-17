#!/usr/bin/env python3
from textwrap import dedent

import hailtop.batch as hb
from os.path import join
import os

from cpg_pipes import resources
from cpg_pipes.jobs import wrap_command

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
BWAMEM2_IMAGE = f'{AR_REPO}/alignment:v4'
PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')


def _index_bwa_job(b: hb.Batch):
    reference = b.read_input_group(**resources.REF_D)

    j = b.new_job('Index BWA')
    j.image(BWAMEM2_IMAGE)
    total_cpu = 32
    j.cpu(total_cpu)
    j.storage('40G')
    j.declare_resource_group(
        bwa_index={e: '{root}.' + e for e in resources.BWAMEM2_INDEX_EXTS}
    )
    j.command(wrap_command(f"""\
    set -o pipefail
    set -ex
    
    bwa-mem2 index {reference.base} -p {j.bwa_index}
    
    df -h; pwd; ls | grep -v proc | xargs du -sh
    """))
    b.write_output(j.bwa_index, resources.BWAMEM2_INDEX_PREFIX)
    return j


def _test_bwa_job(b: hb.Batch):
    bwa_reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
        sa=REF_FASTA + '.sa',
        amb=REF_FASTA + '.amb',
        bwt=REF_FASTA + '.bwt',
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
        o123=REF_FASTA + '.0123',
        bwa2bit64=REF_FASTA + '.bwt.2bit.64',
    )

    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test BWA')
    j.image(BWAMEM2_IMAGE)
    total_cpu = 16
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
    rg_line = f'@RG\\tID|:{sn}\\tSM:~{sn}'
    # bwa index {reference.base}
    use_bazam = True

    sorted_bam = f'$(dirname {j.output_cram.cram})/sorted.bam'
    j.command(
        dedent(
            f"""
    set -o pipefail
    set -ex
    
    (while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &
    
    bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
        -R '{rg_line}' {bwa_reference.base} {fq1} - | \\
    samtools sort -T $(dirname {j.output_cram.cram})/samtools-sort-tmp -Obam -o {sorted_bam}
    
    picard MarkDuplicates I={sorted_bam} O=/dev/stdout M={j.duplicate_metrics} \\
        ASSUME_SORT_ORDER=coordinate | \\
    samtools view -@30 -T {bwa_reference.base} -O cram -o {j.output_cram.cram}
    
    samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}
    
    df -h; pwd; du -sh $(dirname {j.output_cram.cram})
    """
        )
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
b = hb.Batch(backend=backend, name='Create BWA index')
j1 = _index_bwa_job(b)
j2 = _test_bwa_job(b)
j2.depends_on(j1)
b.run(wait=False)
