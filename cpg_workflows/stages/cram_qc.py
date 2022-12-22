"""
Stages that generates and summarises CRAM QC.
"""
import logging
import dataclasses
from typing import Callable, Optional

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import somalier
from cpg_workflows.jobs.multiqc import multiqc
from cpg_workflows.jobs.picard import (
    picard_wgs_metrics,
    picard_collect_metrics,
    picard_hs_metrics,
)
from cpg_workflows.jobs.samtools import samtools_stats
from cpg_workflows.jobs.verifybamid import verifybamid
from cpg_workflows.stages.align import Align
from cpg_workflows.targets import Sample, Dataset
from cpg_workflows.utils import exists
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    DatasetStage,
    StageInputNotFoundError,
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
    optional: bool = False


def qc_functions() -> list[Qc]:
    """
    QC functions and their outputs for MultiQC aggregation
    """
    if get_config()['workflow'].get('skip_qc', False) is True:
        return []

    qcs = [
        Qc(func=somalier.extract, outs={'somalier': None}),
        Qc(
            func=verifybamid,
            outs={'verify_bamid': QcOut('.verify-bamid.selfSM', 'verifybamid/selfsm')},
        ),
        Qc(
            func=samtools_stats,
            outs={'samtools_stats': QcOut('.samtools-stats', 'samtools/stats')},
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


@stage(required_stages=Align)
class CramQC(SampleStage):
    """
    Calling tools that process CRAM for QC purposes.
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        outs = {}
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
        cram_path = inputs.as_path(sample, Align, 'cram')
        crai_path = inputs.as_path(sample, Align, 'crai')

        jobs = []
        for qc in qc_functions():
            out_path_kwargs = {
                f'out_{key}_path': self.expected_outputs(sample)[key]
                for key in qc.outs.keys()
            }
            if qc.func:
                j = qc.func(  # type: ignore
                    get_batch(),
                    CramPath(cram_path, crai_path),
                    job_attrs=self.get_job_attrs(sample),
                    overwrite=sample.forced,
                    **out_path_kwargs,
                )
                if j:
                    jobs.append(j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)


@stage(required_stages=[CramQC], forced=True)
class SomalierPedigree(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}
        h = dataset.alignment_inputs_hash()
        prefix = dataset.prefix() / 'somalier' / 'cram' / h
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': dataset.web_prefix() / 'cram-somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Checks calls job from the pedigree module
        """
        verifybamid_by_sid = {}
        somalier_path_by_sid = {}
        for sample in dataset.get_samples():
            if get_config().get('somalier', {}).get('exclude_high_contamination'):
                verify_bamid_path = inputs.as_path(
                    stage=CramQC, target=sample, key='verify_bamid'
                )
                if not exists(verify_bamid_path):
                    logging.warning(
                        f'VerifyBAMID results {verify_bamid_path} do not exist for '
                        f'{sample}, somalier pedigree estimations might be affected'
                    )
                else:
                    verifybamid_by_sid[sample.id] = verify_bamid_path
            somalier_path = inputs.as_path(stage=CramQC, target=sample, key='somalier')
            somalier_path_by_sid[sample.id] = somalier_path

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        if any(s.pedigree for s in dataset.get_samples()):
            expected_ped_path = dataset.write_ped_file(
                self.expected_outputs(dataset)['expected_ped']
            )
            jobs = somalier.pedigree(
                b=get_batch(),
                dataset=dataset,
                expected_ped_path=expected_ped_path,
                somalier_path_by_sid=somalier_path_by_sid,
                verifybamid_by_sid=verifybamid_by_sid,
                out_samples_path=self.expected_outputs(dataset)['samples'],
                out_pairs_path=self.expected_outputs(dataset)['pairs'],
                out_html_path=html_path,
                out_html_url=html_url,
                out_checks_path=self.expected_outputs(dataset)['checks'],
                job_attrs=self.get_job_attrs(dataset),
                send_to_slack=True,
            )
            return self.make_outputs(
                dataset, data=self.expected_outputs(dataset), jobs=jobs
            )
        else:
            return self.make_outputs(dataset, skipped=True)


@stage(
    required_stages=[
        CramQC,
        SomalierPedigree,
    ],
    forced=True,
)
class CramMultiQC(DatasetStage):
    """
    Run MultiQC to aggregate CRAM QC stats.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}
        h = dataset.alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'cram' / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'cram' / h / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'cram' / h / '.checks',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Call a function from the `jobs` module using inputs from `cramqc`
        and `somalier` stages.
        """
        if get_config()['workflow'].get('skip_qc', False) is True:
            return self.make_outputs(dataset)

        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        checks_path = self.expected_outputs(dataset)['checks']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        paths = []
        try:
            somalier_samples = inputs.as_path(dataset, SomalierPedigree, key='samples')
            somalier_pairs = inputs.as_path(dataset, SomalierPedigree, key='pairs')
        except StageInputNotFoundError:
            pass
        else:
            paths = [
                somalier_samples,
                somalier_pairs,
            ]

        ending_to_trim = set()  # endings to trim to get sample names
        modules_to_trim_endings = set()

        for sample in dataset.get_samples():
            for qc in qc_functions():
                for key, out in qc.outs.items():
                    if not out:
                        continue
                    try:
                        path = inputs.as_path(sample, CramQC, key)
                    except StageInputNotFoundError:  # allow missing inputs
                        logging.warning(
                            f'Output CramQc/"{key}" not found for {sample}, '
                            f'it will be silently excluded from MultiQC'
                        )
                        continue
                    modules_to_trim_endings.add(out.multiqc_key)
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sample.id, ''))

        if not paths:
            logging.warning('No CRAM QC found to aggregate with MultiQC')
            return self.make_outputs(dataset)

        jobs = multiqc(
            get_batch(),
            tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'cram',
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            dataset=dataset,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=checks_path,
            job_attrs=self.get_job_attrs(dataset),
            sample_id_map=dataset.rich_id_map(),
            label='CRAM',
            extra_config={'table_columns_visible': {'FastQC': False}},
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )
