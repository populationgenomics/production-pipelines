"""
Stage that generates a CRAM file.
"""
import logging
import dataclasses
from typing import Callable, Optional

from cpg_utils.config import get_config
from cpg_utils import Path
from cpg_utils.workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
)

from cpg_workflows.jobs import align, somalier
from cpg_workflows.jobs.align import MissingAlignmentInputException
from cpg_workflows.jobs.verifybamid import verifybamid
from cpg_workflows.jobs.picard import (
    picard_wgs_metrics,
    picard_collect_metrics,
    picard_hs_metrics,
)


@dataclasses.dataclass
class QcOut:
    """QC output file"""

    suf: str
    multiqc_key: str


@dataclasses.dataclass
class Qc:
    """QC function definition, and corresponding outputs"""

    func: Optional[Callable]
    outs: dict[str, QcOut | None]


def qc_functions() -> list[Qc]:
    """
    QC functions and their outputs for MultiQC aggregation
    """
    if get_config()['workflow'].get('skip_qc', False) is True:
        return []

    qcs = [
        Qc(
            func=None,
            outs={
                'markduplicates_metrics': QcOut(
                    '.markduplicates-metrics', 'picard/markdups'
                )
            },
        ),
        Qc(func=somalier.extract, outs={'somalier': None}),
        Qc(
            func=verifybamid,
            outs={'verify_bamid': QcOut('.verify-bamid.selfSM', 'verifybamid/selfsm')},
        ),
        Qc(
            func=picard_collect_metrics,
            outs={
                'alignment_summary_metrics': QcOut(
                    '.alignment_summary_metrics', 'picard/alignment_metrics'
                ),
                'base_distribution_by_cycle_metrics': QcOut(
                    '.base_distribution_by_cycle_metrics',
                    'picard/basedistributionbycycle',
                ),
                'insert_size_metrics': QcOut(
                    '.insert_size_metrics', 'picard/insertsize'
                ),
                'quality_by_cycle_metrics': QcOut(
                    '.quality_by_cycle_metrics', 'picard/quality_by_cycle'
                ),
                'quality_yield_metrics': QcOut(
                    '.quality_yield_metrics', 'picard/quality_yield_metrics'
                ),
            },
        ),
    ]

    sequencing_type = get_config()['workflow']['sequencing_type']
    if sequencing_type == 'genome':
        qcs.append(
            Qc(
                func=picard_wgs_metrics,
                outs={
                    'picard_wgs_metrics': QcOut(
                        '.picard-wgs-metrics', 'picard/wgs_metrics'
                    )
                },
            )
        )
    if sequencing_type == 'exome':
        qcs.append(
            Qc(
                func=picard_hs_metrics,
                outs={
                    'picard_hs_metrics': QcOut('.picard-hs-metrics', 'picard/hsmetrics')
                },
            )
        )
    return qcs


@stage(analysis_type='cram')
class Align(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        outs = {
            'cram': sample.make_cram_path().path,
        }
        for qc in qc_functions():
            for key, out in qc.outs.items():
                if key == 'somalier':
                    outs[key] = sample.make_cram_path().somalier_path
                elif out:
                    outs[key] = (
                        sample.dataset.prefix() / 'qc' / key / f'{sample.id}{out.suf}'
                    )
        return outs

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
        """
        jobs = []
        try:
            align_jobs = align.align(
                b=self.b,
                sample=sample,
                output_path=sample.make_cram_path(),
                out_markdup_metrics_path=self.expected_outputs(sample).get(
                    'markduplicates_metrics'
                ),
                job_attrs=self.get_job_attrs(sample),
                overwrite=not get_config()['workflow'].get('check_intermediates'),
            )
        except MissingAlignmentInputException:
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logging.error(f'No alignment inputs, skipping sample {sample}')
                sample.active = False
                return self.make_outputs(sample, skipped=True)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, error_msg=f'No alignment input found'
                )
        else:
            jobs.extend(align_jobs)

        for qc in qc_functions():
            out_path_kwargs = {
                f'out_{key}_path': self.expected_outputs(sample)[key]
                for key in qc.outs.keys()
            }
            if qc.func:
                j = qc.func(  # type: ignore
                    self.b,
                    sample.make_cram_path(),
                    job_attrs=self.get_job_attrs(sample),
                    overwrite=not get_config()['workflow'].get('check_intermediates'),
                    **out_path_kwargs,
                )
                if j:
                    if align_jobs:
                        j.depends_on(*align_jobs)
                    jobs.append(j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
