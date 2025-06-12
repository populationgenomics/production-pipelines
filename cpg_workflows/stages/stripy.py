"""
Stage to run STR analysis with STRipy-pipeline.

See https://gitlab.com/andreassh/stripy-pipeline
"""

from typing import Any

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import stripy
from cpg_workflows.stages.align import Align
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


def _stripy_config_retrieve(key: str, dataset: str | None = None) -> Any:
    """
    Retrieve the STRipy config for a given key. If a dataset is provided, it will first
    check for a dataset-specific config before falling back to the global config.

    **NOT to be confused with the stripy config.json file used by the stripy job**
    """
    if dataset and config_retrieve(['stripy', dataset, key], default=None) is not None:
        return config_retrieve(['stripy', dataset, key], default=None)
    else:
        return config_retrieve(['stripy', key], default=None)


def _get_target_loci(dataset: str) -> str:
    """
    Get the target loci which stripy should run on.
    Uses defaults + dataset-specific additional loci if available.
    """
    target_loci = _stripy_config_retrieve('target_loci')
    if dataset_specific_loci := _stripy_config_retrieve('target_loci', dataset):
        if _stripy_config_retrieve('target_loci_override', dataset) is True:
            target_loci = dataset_specific_loci
        else:
            target_loci = ','.join(
                set(target_loci.split(',')) | set(dataset_specific_loci.split(',')),
            )
    return target_loci


def _update_meta(output_path: str) -> dict[str, Any]:
    """
    Add the detected outlier loci to the analysis meta
    """
    from cloudpathlib import CloudPath

    # Munge html path into log path (As far as I can know I can not pass to
    # output paths to one analysis object?)
    log_path = output_path.replace('-web/', '-analysis/').replace('.html', '.log.txt')

    outlier_loci = {}
    with CloudPath(log_path).open() as f:
        for line in f:
            path, symbol, score = line.strip().split('\t')
            if not score.isdigit():
                continue
            if int(score) > 0:
                outlier_loci[symbol] = score

    return {
        'outlier_loci': outlier_loci,
        'outliers_detected': bool(outlier_loci),
        'log_path': log_path,
    }


@stage(
    required_stages=Align,
    analysis_type='web',
    analysis_keys=[
        'stripy_html',
    ],
    update_analysis_meta=_update_meta,
)
class Stripy(SequencingGroupStage):
    """
    Call stripy to run STR analysis on known pathogenic loci.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        # When testing new loci, its useful to save the report to an alternate path
        output_prefix = config_retrieve(['stripy', 'output_prefix'], default='stripy')
        return {
            'stripy_html': sequencing_group.dataset.web_prefix() / output_prefix / f'{sequencing_group.id}.stripy.html',
            'stripy_json': sequencing_group.dataset.analysis_prefix()
            / output_prefix
            / f'{sequencing_group.id}.stripy.json',
            'stripy_log': sequencing_group.dataset.analysis_prefix()
            / output_prefix
            / f'{sequencing_group.id}.stripy.log.txt',
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path = inputs.as_path(sequencing_group, Align, 'cram')
        crai_path = inputs.as_path(sequencing_group, Align, 'crai')

        jobs = []
        j = stripy.stripy(
            b=get_batch(),
            sequencing_group=sequencing_group,
            cram_path=CramPath(cram_path, crai_path),
            target_loci=_get_target_loci(sequencing_group.dataset.name),
            custom_loci_path=_stripy_config_retrieve('custom_loci_path', sequencing_group.dataset.name),
            analysis_type=_stripy_config_retrieve('analysis_type'),
            stripy_config=_stripy_config_retrieve('config', sequencing_group.dataset.name),
            log_path=self.expected_outputs(sequencing_group)['stripy_log'],
            out_path=self.expected_outputs(sequencing_group)['stripy_html'],
            json_path=self.expected_outputs(sequencing_group)['stripy_json'],
            job_attrs=self.get_job_attrs(sequencing_group),
        )
        jobs.append(j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)
