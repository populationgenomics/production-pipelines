"""
Stage that runs FastQC on alignment inputs.
"""
import dataclasses

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import BamPath, FastqPairs
from cpg_workflows.workflow import (
    Sample,
    Dataset,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    DatasetStage,
)
from cpg_workflows.jobs import fastqc
from cpg_workflows.jobs.multiqc import multiqc


@dataclasses.dataclass
class OneFastqc:
    """
    Inputs and outputs for one FASTQC job.
    """

    suffix: str
    input_path: Path
    out_html: Path
    out_zip: Path


def _collect_fastq_outs(sample) -> list[OneFastqc]:
    """
    Collect input and output paths for FASTQC for all paths in alignment inputs.
    """
    sequencing_type = get_config()['workflow']['sequencing_type']
    if not (alignment_input := sample.alignment_input_by_seq_type.get(sequencing_type)):
        # Only running FASTQC if sequencing inputs are available.
        return []

    if get_config()['workflow'].get('check_inputs', True):
        if not alignment_input.exists():
            return []

    prefix = sample.dataset.prefix() / 'qc' / 'fastqc'

    if isinstance(alignment_input, BamPath):
        return [
            OneFastqc(
                '',
                alignment_input.path,
                prefix / (sample.id + '_fastqc.html'),
                prefix / (sample.id + '_fastqc.zip'),
            )
        ]
    elif isinstance(alignment_input, FastqPairs):
        outs: list[OneFastqc] = []
        for lane_i, pair in enumerate(alignment_input):
            lane_suffix = f'_lane{lane_i + 1}' if len(alignment_input) > 1 else ''
            for pair_id in [0, 1]:
                suffix = f'{lane_suffix}_R{pair_id + 1}'
                outs.append(
                    OneFastqc(
                        suffix,
                        pair[pair_id],
                        prefix / f'{sample.id}{suffix}_fastqc.html',
                        prefix / f'{sample.id}{suffix}_fastqc.zip',
                    )
                )
        return outs
    return []


@stage
class FastQC(SampleStage):
    """
    Run FASTQC on all paths in alignment inputs.
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path] | None:
        """
        Generates one FASTQC HTML report per "sequence" path
        (a FASTQ path, or a BAM path depending on the inputs type).
        """
        outs: dict[str, Path] = {}
        for fq_out in _collect_fastq_outs(sample):
            outs |= {
                f'html{fq_out.suffix}': fq_out.out_html,
                f'zip{fq_out.suffix}': fq_out.out_zip,
            }
        return outs

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        if not (fqc_outs := _collect_fastq_outs(sample)):
            return self.make_outputs(sample, skipped=True)

        jobs = []
        for fqc_out in fqc_outs:
            j = fastqc.fastqc(
                b=get_batch(),
                output_html_path=fqc_out.out_html,
                output_zip_path=fqc_out.out_zip,
                input_path=fqc_out.input_path,
                job_attrs=self.get_job_attrs(sample),
                subsample=False,
            )
            j.name = f'{j.name}{fqc_out.suffix}'
            jobs.append(j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)


@stage(required_stages=FastQC, forced=True)
class FastQCMultiQC(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}
        h = dataset.alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'fastqc' / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'fastqc' / h / 'multiqc_data.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        paths = []  # FASTQC zip outputs to parse with MultiQC
        sample_id_map = {}  # MultiQC would use names of FASTQ files as sample names,
        # we need to collect a map to rename them to proper internal/external IDs
        for sample in dataset.get_samples():
            for fqc_out in _collect_fastq_outs(sample):
                paths.append(fqc_out.out_zip)
                fq_name = fqc_out.input_path.name.removesuffix('.gz').split('.')[0]
                sample_id_map[fq_name] = f'{sample.rich_id}{fqc_out.suffix}'

        jobs = multiqc(
            get_batch(),
            tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'fastqc',
            paths=paths,
            dataset=dataset,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(dataset),
            sample_id_map=sample_id_map,
            label='FASTQC',
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )
