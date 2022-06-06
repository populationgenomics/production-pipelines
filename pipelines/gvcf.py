#!/usr/bin/env python3

"""
Batch pipeline to generate GVCF for all samples CRAMs.
"""

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.genotype_sample import GenotypeSample


pipeline = Pipeline('Sample genotyping')
pipeline.run(stages=[GenotypeSample])
