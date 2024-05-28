#!/usr/bin/env python3
"""
Production version of workflows.
"""
from cpg_workflows.stages.cram_qc import CramMultiQC
from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
from cpg_workflows.stages.large_cohort import AncestryPlots, Frequencies, LoadVqsr
from cpg_workflows.workflow import StageDecorator

from cpg_workflows.stages.aip import CreateAIPHTML, GenerateSeqrFile, ValidateMOI
from cpg_workflows.stages.exomiser import ExomiserSeqrTSV, RunExomiser
from cpg_workflows.stages.fastqc import FastQCMultiQC
from cpg_workflows.stages.fraser import Fraser
from cpg_workflows.stages.gatk_sv.gatk_sv_multisample_1 import (
    FilterBatch,
    GenotypeBatch,
    MergeBatchSites,
)
from cpg_workflows.stages.gatk_sv.gatk_sv_multisample_2 import MtToEsSv
from cpg_workflows.stages.gatk_sv.gatk_sv_single_sample import CreateSampleBatches
from cpg_workflows.stages.gcnv import AnnotateCohortgCNV, AnnotateDatasetCNV, MtToEsCNV
from cpg_workflows.stages.happy_validation import (
    ValidationHappyOnVcf,
    ValidationMtToVcf,
    ValidationParseHappy,
)
from cpg_workflows.stages.mito import MitoReport
from cpg_workflows.stages.outrider import Outrider
from cpg_workflows.stages.seqr_loader import AnnotateDataset, DatasetVCF, MtToEs
from cpg_workflows.stages.stripy import Stripy


WORKFLOWS: dict[str, list[StageDecorator]] = {
    'aip': [ValidateMOI, CreateAIPHTML, GenerateSeqrFile],
    'exomiser': [RunExomiser, ExomiserSeqrTSV],
    'pre_alignment': [FastQCMultiQC],
    'seqr_loader': [
        DatasetVCF,
        AnnotateDataset,
        MtToEs,
        GvcfMultiQC,
        CramMultiQC,
        Stripy,
        MitoReport,
    ],
    'validation': [ValidationMtToVcf, ValidationHappyOnVcf, ValidationParseHappy],
    'large_cohort': [LoadVqsr, Frequencies, AncestryPlots, GvcfMultiQC, CramMultiQC],
    'gatk_sv_singlesample': [CreateSampleBatches],
    'gatk_sv_multisample_1': [FilterBatch, GenotypeBatch],
    'gatk_sv_sandwich': [MergeBatchSites],  # stage to run between FilterBatch & GenotypeBatch
    'gatk_sv_multisample_2': [MtToEsSv],
    'rare_disease_rnaseq': [Outrider, Fraser],
    'gcnv': [AnnotateCohortgCNV, AnnotateDatasetCNV, MtToEsCNV],
}