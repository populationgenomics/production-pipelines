#!/usr/bin/env python3

"""
Batch pipeline to run QC on all CRAMs. Outputs MultiQC reports for each dataset, 
and pedigree checks.
"""

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.multiqc import CramMultiQC


pipeline = Pipeline('CRAM QC')
pipeline.run(
    stages=[CramMultiQC],
    force_all_implicit_stages=True,
)
