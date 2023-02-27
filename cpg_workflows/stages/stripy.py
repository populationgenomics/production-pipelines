"""
Stage to run STR analysis with STRipy-pipeline.

See https://gitlab.com/andreassh/stripy-pipeline
"""
from typing import Any

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import stripy
from cpg_workflows.stages.align import Align
from cpg_workflows.targets import Sample
from cpg_workflows.utils import exists
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)


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
            if int(score) > 0:
                outlier_loci[symbol] = score

    return {
        'outlier_loci': outlier_loci,
        'outliers_detected': bool(outlier_loci),
        'log_path': log_path
    }


@stage(
    required_stages=Align,
    analysis_type='web',
    analysis_keys=['stripy_html', ],
    update_analysis_meta=_update_meta,
)
class Stripy(SampleStage):
    """
    Call stripy to run STR analysis on known pathogenic loci.
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        return {
            'stripy_html': sample.dataset.web_prefix()
            / 'stripy'
            / f'{sample.id}.stripy.html',
            'stripy_json': sample.dataset.analysis_prefix()
            / 'stripy'
            / f'{sample.id}.stripy.json',
            'stripy_log': sample.dataset.analysis_prefix()
            / 'stripy'
            / f'{sample.id}.stripy.log.txt',
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        cram_path = inputs.as_path(sample, Align, 'cram')
        crai_path = inputs.as_path(sample, Align, 'crai')

        jobs = []
        j = stripy.stripy(
            b=get_batch(),
            sample=sample,
            cram_path=CramPath(cram_path, crai_path),
            target_loci=get_config()['stripy']['target_loci'],
            log_path=self.expected_outputs(sample)['stripy_log'],
            analysis_type=get_config()['stripy']['analysis_type'],
            out_path=self.expected_outputs(sample)['stripy_html'],
            json_path=self.expected_outputs(sample)['stripy_json'],
            job_attrs=self.get_job_attrs(sample),
        )
        jobs.append(j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
