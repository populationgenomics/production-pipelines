#!/usr/bin/env python3

"""
Benchmarking DRAGMAP alignment and GATK-DRAGEN variant calling.
"""

from os.path import join
import click
import logging

from cpg_pipes import benchmark
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.jobs.align import Aligner, align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


DATASET = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():
    pipe = Pipeline(
        analysis_dataset='fewgenomes',
        name='benchmark_dragen',
        output_version='v0',
        namespace=NAMESPACE,
        description='Benchmark DRAGMAP: full samples',
        check_smdb_seq=False,
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }
    for sample_name, inp in inputs.items():
        cram_path = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}.cram'
        align_j = align(
            pipe.b,
            alignment_input=inp,
            output_path=cram_path,
            sample_name=sample_name,
            dataset_name='benchmark',
            aligner=Aligner.DRAGMAP,
        )
        produce_gvcf(
            pipe.b,
            output_path=f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}.g.vcf.gz',
            sample_name=sample_name,
            dataset_name='benchmark',
            cram_path=cram_path,
            number_of_intervals=10,
            tmp_bucket=join(benchmark.BENCHMARK_BUCKET, 'tmp'),
            depends_on=[align_j],
            smdb=pipe.get_db(),
            dragen_mode=True,
        )

    pipe.submit_batch()


if __name__ == '__main__':
    main()
