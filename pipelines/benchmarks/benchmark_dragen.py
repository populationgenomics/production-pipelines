#!/usr/bin/env python3
from os.path import join
import click
import logging

from cpg_pipes import benchmark
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.pipeline import Pipeline
from cpg_pipes.jobs.align import Aligner, MarkDupTool, \
    align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


PROJECT = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='benchmark_dragen',
        output_version='v0',
        namespace=NAMESPACE,
        title='Benchmark DRAGMAP: full samples',
        check_smdb_seq_existence=False,
    )
    tiny_inputs = {
        'TINY_FQ': benchmark.tiny_fq,
        'TINY_CRAM': benchmark.tiny_cram,
    }
    for sample_name, inp in tiny_inputs.items():
        cram_path = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}.cram'
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
            output_path=f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}.g.vcf.gz',
            sample_name=sample_name,
            project_name='benchmark',
            cram_path=cram_path,
            number_of_intervals=50,
            tmp_bucket=join(benchmark.BENCHMARK_BUCKET, 'tmp'),
            depends_on=[align_j],
            smdb=pipe.db,
            dragen_mode=True,
        )
    
    pipe.run()


if __name__ == '__main__':
    main()
