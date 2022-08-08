#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.somalier import CramSomalierAncestry, CramSomalierPedigree


pipeline = Pipeline('Somalier')
pipeline.run(
    stages=[
        CramSomalierPedigree,
    ],
)
