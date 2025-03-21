#!/usr/bin/env python3

"""
Entry point to run workflows.
"""

import os

import click
import coloredlogs

from cpg_utils import to_path
from cpg_utils.config import set_config_paths
from cpg_workflows import defaults_config_path
from cpg_workflows.stages.clinvarbitration import PackageForRelease
from cpg_workflows.stages.cram_qc import CramMultiQC
from cpg_workflows.stages.exomiser import ExomiserSeqrTSV, ExomiserVariantsTSV, RegisterSingleSampleExomiserResults
from cpg_workflows.stages.fastqc import FastQCMultiQC
from cpg_workflows.stages.fraser import Fraser
from cpg_workflows.stages.gatk_sv.gatk_sv_multisample import (
    FilterBatch,
    GenotypeBatch,
    MtToEsSv,
    SplitAnnotatedSvVcfByDataset,
)
from cpg_workflows.stages.gatk_sv.gatk_sv_single_sample import CreateSampleBatches
from cpg_workflows.stages.gcnv import AnnotateCohortgCNV, AnnotateDatasetCNV, MtToEsCNV, SplitAnnotatedCnvVcfByDataset
from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
from cpg_workflows.stages.happy_validation import ValidationHappyOnVcf, ValidationMtToVcf
from cpg_workflows.stages.large_cohort import AncestryPlots, Frequencies, LoadVqsr
from cpg_workflows.stages.mito import MitoReport
from cpg_workflows.stages.outrider import Outrider
from cpg_workflows.stages.rd_combiner import (
    AnnotateCohortSmallVariantsWithHailQuery,
    AnnotateDatasetSmallVariantsWithHailQuery,
    AnnotateFragmentedVcfWithVep,
    ConcatenateVcfFragmentsWithGcloud,
    CreateDenseMtFromVdsWithHail,
    CreateVdsFromGvcfsWithHailCombiner,
    ExportMtAsEsIndex,
    GatherTrainedVqsrSnpTranches,
    RunTrainedIndelVqsrOnCombinedVcf,
    SubsetMatrixTableToDatasetUsingHailQuery,
    TrainVqsrSnpTranches,
)
from cpg_workflows.stages.seqr_loader import AnnotateDataset, DatasetVCF, MtToEs
from cpg_workflows.stages.seqr_loader_long_read.bam_to_cram import ConvertPacBioBamToCram
from cpg_workflows.stages.seqr_loader_long_read.long_read_snps_indels_annotation import MtToEsLrSNPsIndels
from cpg_workflows.stages.seqr_loader_long_read.long_read_sv_annotation import MtToEsLrSv
from cpg_workflows.stages.stripy import Stripy
from cpg_workflows.stages.talos import MakePhenopackets, MinimiseOutputForSeqr, UploadTalosHtml, ValidateMOI
from cpg_workflows.stages.talos_prep.talos_prep import SquashMtIntoTarball
from cpg_workflows.workflow import StageDecorator, run_workflow

WORKFLOWS: dict[str, list[StageDecorator]] = {
    'clinvarbitration': [PackageForRelease],
    'talos': [MakePhenopackets, ValidateMOI, UploadTalosHtml, MinimiseOutputForSeqr],
    'talos_prep': [SquashMtIntoTarball],
    'exomiser': [ExomiserSeqrTSV, ExomiserVariantsTSV, RegisterSingleSampleExomiserResults],
    'long_read_snps_indels_annotation': [MtToEsLrSNPsIndels],
    'long_read_sv_annotation': [MtToEsLrSv],
    'pre_alignment': [FastQCMultiQC],
    'rd_combiner': [
        CreateVdsFromGvcfsWithHailCombiner,
        CreateDenseMtFromVdsWithHail,
        ConcatenateVcfFragmentsWithGcloud,
        GatherTrainedVqsrSnpTranches,
        TrainVqsrSnpTranches,
        RunTrainedIndelVqsrOnCombinedVcf,
        AnnotateFragmentedVcfWithVep,
        AnnotateCohortSmallVariantsWithHailQuery,
        SubsetMatrixTableToDatasetUsingHailQuery,
        AnnotateDatasetSmallVariantsWithHailQuery,
        ExportMtAsEsIndex,
    ],
    'seqr_loader': [
        DatasetVCF,
        AnnotateDataset,
        MtToEs,
        GvcfMultiQC,
        CramMultiQC,
        Stripy,
        MitoReport,
    ],
    'seqr_loader_long_read': [
        ConvertPacBioBamToCram,
    ],
    'validation': [ValidationMtToVcf, ValidationHappyOnVcf],
    'large_cohort': [LoadVqsr, Frequencies, AncestryPlots, GvcfMultiQC, CramMultiQC],
    'gatk_sv_singlesample': [CreateSampleBatches],
    'gatk_sv_multisample': [FilterBatch, GenotypeBatch, MtToEsSv, SplitAnnotatedSvVcfByDataset],
    'rare_disease_rnaseq': [Outrider, Fraser],
    'gcnv': [AnnotateCohortgCNV, AnnotateDatasetCNV, MtToEsCNV, SplitAnnotatedCnvVcfByDataset],
}


@click.command(no_args_is_help=True)
@click.argument('workflow', required=False)
@click.option(
    '--config',
    'config_paths',
    multiple=True,
    help='Add configuration files to the files specified $CPG_CONFIG_PATH.'
    'Configs are merged left to right, meaning the rightmost file has the'
    'highest priority.',
)
@click.option(
    '--list-workflows',
    'list_workflows',
    is_flag=True,
    help='Only list possible values for WORKFLOW (and available last stages)',
)
@click.option(
    '--list-last-stages',
    'list_last_stages',
    is_flag=True,
    help='Only list possible end stages for a workflow, that can be specified with `workflow/last_stages` in config',
)
@click.option(
    '--dry-run',
    'dry_run',
    is_flag=True,
    help='Dry run: do not actually communicate with Metamist or Hail Batch, '
    'instead only print a final config and stages to be run',
)
@click.option(
    '--verbose',
    'verbose',
    is_flag=True,
)
def main(
    workflow: str,
    config_paths: list[str],
    list_workflows: bool,
    list_last_stages: bool,
    dry_run: bool,
    verbose: bool,
):
    """
    Run a Hail Batch workflow specified as a positional command line argument [WORKFLOW]
    """
    fmt = '%(asctime)s %(levelname)s (%(pathname)s %(lineno)s): %(message)s'
    coloredlogs.install(level='DEBUG' if verbose else 'INFO', fmt=fmt)

    if not workflow and not list_workflows:
        click.echo('You must specify WORKFLOW as a first positional command line argument.')
    if not workflow or list_workflows or workflow == 'list':
        click.echo('Available values for WORKFLOW (and corresponding last stages):')
        for wfl, last_stages in WORKFLOWS.items():
            click.echo(f'\t{wfl} ({", ".join(s.__name__ for s in last_stages)})')
        return

    if list_last_stages:
        click.echo(
            f'Available last stages that can be listed with '
            f'workflow/last_stages for the current workflow "{workflow}":',
        )
        click.echo(f'{", ".join(s.__name__ for s in WORKFLOWS[workflow])}')
        return

    wfl_conf_path = to_path(__file__).parent / f'configs/defaults/{workflow}.toml'
    assert wfl_conf_path.exists(), wfl_conf_path

    for path in config_paths:
        assert to_path(path).exists(), path

    config_paths = os.environ['CPG_CONFIG_PATH'].split(',') + list(config_paths)
    # Assuming the defaults is already loaded in __init__.py:
    assert to_path(config_paths[0]) == defaults_config_path
    # Inserting after the "defaults" config, but before user configs:
    set_config_paths(config_paths[:1] + [str(wfl_conf_path)] + config_paths[1:])

    run_workflow(stages=WORKFLOWS[workflow], dry_run=dry_run)


if __name__ == '__main__':
    main()
