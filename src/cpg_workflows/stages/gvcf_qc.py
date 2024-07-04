"""
Stages that generates and summarises GVCF QC.
"""

import logging
from typing import Any

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs.happy import happy
from cpg_workflows.jobs.multiqc import multiqc
from cpg_workflows.jobs.picard import vcf_qc
from cpg_workflows.stages.genotype import Genotype
from cpg_workflows.targets import Dataset
from cpg_workflows.workflow import (
    DatasetStage,
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageInputNotFoundError,
    StageOutput,
    stage,
)


@stage(required_stages=Genotype)
class GvcfQC(SequencingGroupStage):
    """
    Calling tools that process GVCF for QC purposes.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate a GVCF and corresponding TBI index, as well as QC.
        """
        outs: dict[str, Path] = {}
        if not config_retrieve(['workflow', 'skip_qc'], False):
            qc_prefix = sequencing_group.dataset.prefix() / 'qc' / sequencing_group.id
            outs |= {
                'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
                'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
            }
        return outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        gvcf_path = inputs.as_path(sequencing_group, Genotype, 'gvcf')

        j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(get_batch()),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_summary_path=self.expected_outputs(sequencing_group)['qc_summary'],
            output_detail_path=self.expected_outputs(sequencing_group)['qc_detail'],
            overwrite=sequencing_group.forced,
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=[j])


@stage(required_stages=Genotype)
class GvcfHappy(SequencingGroupStage):
    """
    Run Happy to validate a GVCF for samples where a truth callset is available.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> Path | None:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        if sequencing_group.participant_id not in config_retrieve(['validation', 'sample_map'], {}):
            return None
        return sequencing_group.dataset.prefix() / 'qc' / 'gvcf' / 'hap.py' / f'{sequencing_group.id}.summary.csv'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        gvcf_path = inputs.as_path(sequencing_group, Genotype, 'gvcf')

        jobs = happy(
            b=get_batch(),
            sequencing_group=sequencing_group,
            vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(get_batch()),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_path=self.expected_outputs(sequencing_group),
        )

        if not jobs:
            return self.make_outputs(sequencing_group)
        else:
            return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs)


def _update_meta(output_path: str) -> dict[str, Any]:
    import json

    from cloudpathlib import CloudPath

    with CloudPath(output_path).open() as f:
        d = json.load(f)
    return {'multiqc': d['report_general_stats_data']}


@stage(
    required_stages=[
        GvcfQC,
        GvcfHappy,
    ],
    analysis_type='qc',
    analysis_keys=['json'],
    update_analysis_meta=_update_meta,
)
class GvcfMultiQC(DatasetStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        # get the unique hash for these Sequencing Groups
        sg_hash = dataset.alignment_inputs_hash()
        return {
            'html': dataset.web_prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc_data.json',
            'checks': dataset.prefix() / 'qc' / 'gvcf' / sg_hash / '.checks',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
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

        for sequencing_group in dataset.get_sequencing_groups():
            for _stage, key in [
                (GvcfQC, 'qc_detail'),
                (GvcfHappy, None),
            ]:
                try:
                    path = inputs.as_path(sequencing_group, _stage, key)
                except StageInputNotFoundError:  # allow missing inputs
                    if _stage != GvcfHappy:
                        logging.warning(
                            f'Output {_stage.__name__}/"{key}" not found for {sequencing_group}, '
                            f'it will be silently excluded from MultiQC',
                        )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sequencing_group.id, ''))

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
            sequencing_group_id_map=dataset.rich_id_map(),
            extra_config={'table_columns_visible': {'Picard': True}},
            label='GVCF',
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
