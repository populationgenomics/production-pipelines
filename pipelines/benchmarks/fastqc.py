#!/usr/bin/env python3

"""
Batch pipeline to run FastQC.
"""

import logging
import click
from cloudpathlib import CloudPath

from cpg_pipes import benchmark
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.jobs.fastqc import fastqc

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_DATASET = 'fewgenomes'
NAMESPACE = 'main'


@click.command()
def main():  # pylint: disable=missing-function-docstring
    pipe = Pipeline(
        analysis_dataset='fewgenomes',
        name='run_qc',
        output_version='v0',
        namespace=NAMESPACE,
        description='Run QC',
        check_smdb_seq=False,
    )
    inputs = {
        'NA12878': benchmark.na12878fq,
    }

    for sample_name, inp in inputs.items():
        fastqc(
            pipe.b,
            output_fpath=(
                CloudPath(benchmark.BENCHMARK_BUCKET) / 
                'outputs' /
                'fastqc' /
                f'{sample_name}.html'
            ),
            alignment_input=inp,
            sample_name=sample_name,
            dataset_name='FastQC',
        )
    
    pipe.submit_batch()


if __name__ == '__main__':
    main()
