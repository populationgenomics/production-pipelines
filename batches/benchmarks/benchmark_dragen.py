#!/usr/bin/env python3
from os.path import join

import click
import logging

from cpg_production_pipelines.jobs.haplotype_caller import produce_gvcf
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.align import AlignmentInput, Aligner, MarkDupTool, \
    align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_PROJECT = 'fewgenomes'
NAMESPACE = 'main'
BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test/benchmark'
TOY_INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs/test'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='benchmark_dragen',
        output_version='v0',
        namespace=NAMESPACE,
        title='Benchmark DRAGMAP - full samples',
        smdb_check_existence=False,
    )
    # This samples crashes with BWA-MEM2 after 30min
    # cpg13326 = AlignmentInput(
    #     fqs1=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r1.fastq.gz'],
    #     fqs2=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r2.fastq.gz'],
    # )
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

    for sample_name, inp in giab_inputs.items():
    # for sample_name, inp in tiny_inputs.items():
        cram_path = f'{BENCHMARK_BUCKET}/{sample_name}.cram'
        align_j = align(
            pipe.b,
            alignment_input=inp,
            output_path=cram_path,
            sample_name=sample_name,
            project_name='benchmark',
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.BIOBAMBAM,
            extra_label='biobambam',
        )

        produce_gvcf(
            pipe.b,
            output_path=f'{BENCHMARK_BUCKET}/{sample_name}.g.vcf.gz',
            sample_name=sample_name,
            project_name='benchmark',
            cram_path=cram_path,
            number_of_intervals=50,
            tmp_bucket=join(BENCHMARK_BUCKET, 'tmp'),
            depends_on=[align_j],
            smdb=pipe.db,
            dragen_mode=True,
        )
    
    pipe.run()


if __name__ == '__main__':
    main()
