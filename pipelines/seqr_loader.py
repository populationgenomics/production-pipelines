#!/usr/bin/env python3

"""
Seqr loading pipeline: FASTQ -> ElasticSearch index.
"""

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.seqr_loader import LoadToEs


pipeline = Pipeline(name='Seqr Loader')
pipeline.run(stages=[LoadToEs], force_all_implicit_stages=True)
