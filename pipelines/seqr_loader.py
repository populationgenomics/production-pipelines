#!/usr/bin/env python3

"""
Seqr loading pipeline: FASTQ -> ElasticSearch index.
"""

from cpg_utils.config import get_config
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.seqr_loader import LoadToEs

SUPPORTED_SEQUENCING_TYPES = ['genome']


if seq_type := get_config()['workflow'].get('sequencing_type'):
    if seq_type not in SUPPORTED_SEQUENCING_TYPES:
        raise ValueError(
            f'Unsupported sequencing data type {seq_type.value}. '
            f'Supported types: {[st for st in SUPPORTED_SEQUENCING_TYPES]} '
        )
pipeline = Pipeline(name='Seqr Loader')
pipeline.run(stages=[LoadToEs], force_all_implicit_stages=True)
