#!/usr/bin/env python3

"""
Benchmarking DRAGMAP alignment and GATK-DRAGEN variant calling.
"""

from os.path import join
import click
import logging

from cpg_pipes import benchmark, Namespace
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.jobs.align import Aligner, align
from cpg_pipes.types import SequencingType

logger = logging.getLogger(__file__)


DATASET = 'fewgenomes'
NAMESPACE = Namespace.MAIN


@click.command()
def main():
    pipe = create_pipeline(
        analysis_dataset='fewgenomes',
        name='benchmark_dragen',
        version='v0',
        namespace=NAMESPACE,
        description='Benchmark DRAGMAP: full samples',
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }
    for sample_name, inp in inputs.items():
        cram_path = benchmark.BENCHMARK_BUCKET / 'outputs' / f'{sample_name}.cram'
        align_jobs = align(
            pipe.b,
            alignment_input=inp,
            output_path=cram_path,
            sample_name=sample_name,
            aligner=Aligner.DRAGMAP,
            refs=pipe.refs,
        )
        hc_jobs = produce_gvcf(
            pipe.b,
            output_path=benchmark.BENCHMARK_BUCKET / 'outputs'
            f'{sample_name}.g.vcf.gz',
            sample_name=sample_name,
            cram_path=cram_path,
            scatter_count=10,
            tmp_bucket=benchmark.BENCHMARK_BUCKET / 'tmp',
            dragen_mode=True,
            refs=pipe.refs,
            sequencing_type=SequencingType.GENOME,
        )
        for j in hc_jobs:
            j.depends_on(*align_jobs)

    pipe.run()


if __name__ == '__main__':
    main()
