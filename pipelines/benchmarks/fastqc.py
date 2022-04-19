#!/usr/bin/env python3

"""
Batch pipeline to run FastQC.
"""

import logging
import click

from cpg_pipes import to_path, Namespace
from cpg_pipes import benchmark
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.jobs.fastqc import fastqc

logger = logging.getLogger(__file__)


INPUT_DATASET = 'fewgenomes'
NAMESPACE = Namespace.MAIN


@click.command()
def main():  # pylint: disable=missing-function-docstring
    """
    Entry point
    """
    pipe = create_pipeline(
        analysis_dataset='fewgenomes',
        name='run_qc',
        namespace=NAMESPACE,
        description='Run QC',
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }

    for sample_name, inp in inputs.items():
        prefix = to_path(benchmark.BENCHMARK_BUCKET) / 'outputs' / 'fastqc'
        fastqc(
            pipe.b,
            output_html_path=prefix / f'{sample_name}.html',
            output_zip_path=prefix / f'{sample_name}.zip',
            alignment_input=inp,
            refs=pipe.refs,
        )

    pipe.run()


if __name__ == '__main__':
    main()
