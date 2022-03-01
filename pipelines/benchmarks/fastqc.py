#!/usr/bin/env python3

"""
Batch pipeline to run fastqc
"""

import logging
import click

from cpg_pipes import benchmark
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.jobs.fastqc import fastqc

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_PROJECT = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():  # pylint: disable=missing-function-docstring
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='run_qc',
        output_version='v0',
        namespace=NAMESPACE,
        title='Run QC',
        check_smdb_seq=False,
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
    
    pipe.submit_batch()


if __name__ == '__main__':
    main()
