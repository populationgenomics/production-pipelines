"""
Stages that generates and summarises GVCF QC.
"""
import logging

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs.multiqc import multiqc

from cpg_workflows.stages.genotype import Genotype
from cpg_workflows.targets import Dataset
from cpg_workflows.workflow import (
    Sample,
    stage,
    StageInput,
    StageOutput,
    SampleStage,
    DatasetStage,
    StageInputNotFoundError,
)
from cpg_workflows.jobs.happy import happy
from cpg_workflows.jobs.picard import vcf_qc


@stage(required_stages=Genotype)
class GvcfQC(SampleStage):
    """
    Calling tools that process GVCF for QC purposes.
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Generate a GVCF and corresponding TBI index, as well as QC.
        """
        outs: dict[str, Path] = {}
        if get_config()['workflow'].get('skip_qc', False) is False:
            qc_prefix = sample.dataset.prefix() / 'qc' / sample.id
            outs |= {
                'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
                'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
            }
        return outs

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        gvcf_path = inputs.as_path(sample, Genotype, 'gvcf')

        j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(get_batch()),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sample),
            output_summary_path=self.expected_outputs(sample)['qc_summary'],
            output_detail_path=self.expected_outputs(sample)['qc_detail'],
            overwrite=sample.forced,
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=[j])


@stage(required_stages=Genotype)
class GvcfHappy(SampleStage):
    """
    Run Happy to validate a GVCF for samples where a truth callset is available.
    """

    def expected_outputs(self, sample: Sample) -> Path | None:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        if sample.participant_id not in get_config().get('validation', {}).get(
            'sample_map', {}
        ):
            return None
        return (
            sample.dataset.prefix()
            / 'qc'
            / 'gvcf'
            / 'hap.py'
            / f'{sample.id}.summary.csv'
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        gvcf_path = inputs.as_path(sample, Genotype, 'gvcf')

        jobs = happy(
            b=get_batch(),
            sample=sample,
            vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(get_batch()),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sample),
            output_path=self.expected_outputs(sample),
        )

        if not jobs:
            return self.make_outputs(sample)
        else:
            return self.make_outputs(sample, self.expected_outputs(sample), jobs)


@stage(
    required_stages=[
        GvcfQC,
        GvcfHappy,
    ],
    forced=True,
)
class GvcfMultiQC(DatasetStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}

        h = dataset.alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'gvcf' / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'gvcf' / h / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'gvcf' / h / '.checks',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
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
        ending_to_trim = set()  # endings to trim to get sample names

        for sample in dataset.get_samples():
            for _stage, key in [
                (GvcfQC, 'qc_detail'),
                (GvcfHappy, None),
            ]:
                try:
                    path = inputs.as_path(sample, _stage, key)
                except StageInputNotFoundError:  # allow missing inputs
                    if _stage != GvcfHappy:
                        logging.warning(
                            f'Output {_stage.__name__}/"{key}" not found for {sample}, '
                            f'it will be silently excluded from MultiQC'
                        )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sample.id, ''))

        if not paths:
            logging.warning('No GVCF QC found to aggregate with MultiQC')
            return self.make_outputs(dataset)

        modules_to_trim_endings = {
            'picard/variant_calling_metrics',
            'happy',
        }

        jobs = multiqc(
            get_batch(),
            tmp_prefix=dataset.tmp_prefix() / 'multiqc' / 'gvcf',
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
            extra_config={'table_columns_visible': {'Picard': True}},
            label='GVCF',
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )
