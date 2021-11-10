#!/usr/bin/env python3
import click
import logging

from cpg_production_pipelines import benchmark
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.fastqc import fastqc

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_PROJECT = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='run_qc',
        output_version='v0',
        namespace=NAMESPACE,
        title='Run QC',
        smdb_check_existence=False,
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }

    for sample_name, inp in inputs.items():
        fastqc(
            pipe.b,
            results_bucket=f'{benchmark.BENCHMARK_BUCKET}/outputs/fastqc/{sample_name}',
            alignment_input=inp,
            sample_name=sample_name,
            project_name='FastQC',
        )
    
    pipe.run()


if __name__ == '__main__':
    main()
