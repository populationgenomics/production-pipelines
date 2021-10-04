#!/usr/bin/env python3

import hailtop.batch as hb
import os
import logging
from cpg_production_pipelines.jobs import align
from cpg_production_pipelines.jobs.align import AlignmentInput, Aligner, MarkDupTool

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test/benchmark'
TOY_INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs/test'


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'fewgenomes'
hail_bucket = os.environ.get('HAIL_BUCKET', 'cpg-fewgenomes-test-tmp')

print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='Benchmark alignment')


# This samples crashes with BWA-MEM2 after 30min
cpg13326 = AlignmentInput(
    fqs1=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r1.fastq.gz'],
    fqs2=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r2.fastq.gz'],
)
# This set is 50MB each
tiny_fq = AlignmentInput(
    fqs1=[f'{TOY_INPUTS_BUCKET}/2-699835.L001.R1.n40000.fastq.gz'],
    fqs2=[f'{TOY_INPUTS_BUCKET}/2-699835.L001.R2.n40000.fastq.gz'],
)
tiny_cram_path = f'{TOY_INPUTS_BUCKET}/NA12878-chr21-tiny.cram'
tiny_cram = AlignmentInput(
    bam_or_cram_path=tiny_cram_path,
    index_path=tiny_cram_path + '.crai',
)
tiny_inputs = {
    'TINY_FQ': tiny_fq,
    'TINY_CRAM': tiny_cram,
}

giab_path = 'gs://cpg-reference/validation/giab/cram'
giab_inputs = {
    sn: AlignmentInput(
        bam_or_cram_path=f'{giab_path}/{sn}.cram',
        index_path=f'{giab_path}/{sn}.cram.crai',
    )
    for sn in ['NA12878', 'NA12891', 'NA12892']
}


for sn, inp in giab_inputs.items():
    logger.info(f'Submitting DRAGMAP for {inp}')

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-picard-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.DRAGMAP,
        markdup_tool=MarkDupTool.PICARD,
        extra_label='picard',
    )

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-picard-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.BWA,
        markdup_tool=MarkDupTool.PICARD,
        extra_label='picard',
    )

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-picard-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.BWAMEM2,
        markdup_tool=MarkDupTool.PICARD,
        extra_label='picard',
    )

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-biobambam-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.DRAGMAP,
        markdup_tool=MarkDupTool.BIOBAMBAM,
        extra_label='biobambam',
    )

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-biobambam-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.BWA,
        markdup_tool=MarkDupTool.BIOBAMBAM,
        extra_label='biobambam',
    )

    align.align(
        b=b,
        alignment_input=inp,
        output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-biobambam-from-fq.bam',
        sample_name=sn,
        project_name='Benchmark',
        aligner=Aligner.BWAMEM2,
        markdup_tool=MarkDupTool.BIOBAMBAM,
        extra_label='biobambam',
    )

b.run(wait=False)
