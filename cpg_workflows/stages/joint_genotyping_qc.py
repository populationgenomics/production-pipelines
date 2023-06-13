"""
Stage that summarises QC.
"""
from typing import Any

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    StageInputNotFoundError,
    Cohort,
    SequencingGroupStage,
    get_workflow,
)

from cpg_workflows.jobs.multiqc import multiqc
from cpg_workflows.jobs.picard import vcf_qc
from .joint_genotyping import JointGenotyping
from .. import get_cohort, get_batch
from ..jobs.happy import happy
from ..targets import SequencingGroup


@stage(required_stages=JointGenotyping)
class JointVcfQC(CohortStage):
    """
    QC joint VCF
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Generate a pVCF and a site-only VCF.
        """
        qc_prefix = (
            cohort.analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / get_workflow().output_version
            / 'picard'
        )
        d = {
            'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
            'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
        }
        return d

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        vcf_path = inputs.as_path(target=cohort, stage=JointGenotyping, key='vcf')
        j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=get_batch().read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
                }
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(cohort),
            output_summary_path=self.expected_outputs(cohort)['qc_summary'],
            output_detail_path=self.expected_outputs(cohort)['qc_detail'],
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=JointGenotyping)
class JointVcfHappy(SequencingGroupStage):
    """
    Run Happy to validate validation samples in joint VCF
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> Path | None:
        """
        Parsed by MultiQC: '*.summary.csv'
        https://multiqc.info/docs/#hap.py
        """
        if sequencing_group.participant_id not in get_config().get(
            'validation', {}
        ).get('sample_map', {}):
            return None

        return (
            get_cohort().analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / 'hap.py'
            / f'{get_workflow().output_version}-{sequencing_group.id}.summary.csv'
        )

    def queue_jobs(
        self, sequencing_group: SequencingGroup, inputs: StageInput
    ) -> StageOutput | None:
        """Queue jobs"""
        assert sequencing_group.dataset.cohort
        vcf_path = inputs.as_path(
            target=sequencing_group.dataset.cohort, stage=JointGenotyping, key='vcf'
        )

        jobs = happy(
            b=get_batch(),
            sequencing_group=sequencing_group,
            vcf_or_gvcf=get_batch().read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
                }
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_path=self.expected_outputs(sequencing_group),
        )
        if not jobs:
            return self.make_outputs(sequencing_group)
        else:
            return self.make_outputs(
                sequencing_group, self.expected_outputs(sequencing_group), jobs
            )


def _update_meta(output_path: str) -> dict[str, Any]:
    from cloudpathlib import CloudPath
    import json

    with CloudPath(output_path).open() as f:
        d = json.load(f)
    return {'multiqc': d['report_general_stats_data']}


@stage(
    required_stages=[
        JointVcfQC,
        JointVcfHappy,
    ],
    analysis_type='qc',
    analysis_keys=['json'],
    update_analysis_meta=_update_meta,
)
class JointVcfMultiQC(CohortStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}

        return {
            'html': cohort.analysis_dataset.web_prefix() / 'qc' / 'jc' / 'multiqc.html',
            'json': cohort.analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / get_workflow().output_version
            / 'multiqc_data.json',
            'checks': cohort.analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / get_workflow().output_version
            / '.checks',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        json_path = self.expected_outputs(cohort)['json']
        html_path = self.expected_outputs(cohort)['html']
        checks_path = self.expected_outputs(cohort)['checks']
        if base_url := cohort.analysis_dataset.web_url():
            html_url = str(html_path).replace(
                str(cohort.analysis_dataset.web_prefix()), base_url
            )
        else:
            html_url = None

        paths = []
        ending_to_trim = set()  # endings to trim to get sample names

        paths.append(inputs.as_path(cohort, JointVcfQC, 'qc_detail'))

        for sequencing_group in cohort.get_sequencing_groups():
            try:
                path = inputs.as_path(sequencing_group, JointVcfHappy)
            except StageInputNotFoundError:
                pass
            else:
                paths.append(path)
                ending_to_trim.add(path.name.replace(sequencing_group.id, ''))

        jobs = multiqc(
            get_batch(),
            tmp_prefix=self.tmp_prefix,
            paths=paths,
            ending_to_trim=ending_to_trim,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=checks_path,
            job_attrs=self.get_job_attrs(cohort),
            sequencing_group_id_map=cohort.rich_id_map(),
            extra_config={'table_columns_visible': {'Picard': True}},
            dataset=cohort.analysis_dataset,
            label='Joint VCF',
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
