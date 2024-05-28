#!/usr/bin/env python3
"""
Dev version of workflows.
"""
from cpg_workflows.stages.cram_qc import CramMultiQC
from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
from cpg_workflows.stages.large_cohort import AncestryPlots, Frequencies, LoadVqsr
from cpg_workflows.workflow import StageDecorator

WORKFLOWS: dict[str, list[StageDecorator]] = {
    'large_cohort': [LoadVqsr, Frequencies, AncestryPlots, GvcfMultiQC, CramMultiQC],
}