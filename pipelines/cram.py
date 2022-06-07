#!/usr/bin/env python3

"""
Pipeline to align or realign available sequencing data to produce CRAMs
"""

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.align import Align


pipeline = Pipeline('Alignment')
pipeline.run(stages=[Align])
