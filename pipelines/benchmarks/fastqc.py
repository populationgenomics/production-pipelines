#!/usr/bin/env python3

"""
Batch pipeline to run FastQC.
"""

import logging
import click
from cpg_pipes.storage import to_path
from cpg_pipes import benchmark
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.jobs.fastqc import fastqc

logger = logging.getLogger(__file__)


INPUT_DATASET = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():  # pylint: disable=missing-function-docstring
    """
    Entry point
    """
    pipe = Pipeline(
        analysis_dataset='fewgenomes',
        name='run_qc',
        namespace=NAMESPACE,
        description='Run QC',
        check_smdb_seq=False,
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }

    for sample_name, inp in inputs.items():
        prefix = (
            to_path(benchmark.BENCHMARK_BUCKET) / 
            'outputs' /
            'fastqc'
        )
        fastqc(
            pipe.b,
            output_html_path=prefix / f'{sample_name}.html',
            output_zip_path=prefix / f'{sample_name}.zip',
            alignment_input=inp,
            sample_name=sample_name,
            dataset_name='FastQC',
        )
    
    pipe.submit_batch()


if __name__ == '__main__':
    main()
