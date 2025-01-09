"""
Stage that summarises QC.
"""

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.happy import happy
from cpg_workflows.jobs.multiqc import multiqc
from cpg_workflows.jobs.picard import vcf_qc
from cpg_workflows.stages.joint_genotyping import JointGenotyping
from cpg_workflows.targets import MultiCohort, SequencingGroup
from cpg_workflows.workflow import (
    MultiCohortStage,
    SequencingGroupStage,
    StageInput,
    StageInputNotFoundError,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)


@stage(required_stages=JointGenotyping)
class JointVcfQC(MultiCohortStage):
    """
    QC joint VCF
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        Generate a pVCF and a site-only VCF.
        """
        return {
            'qc_summary': self.prefix / 'picard.variant_calling_summary_metrics',
            'qc_detail': self.prefix / 'picard.variant_calling_detail_metrics',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """

        outputs = self.expected_outputs(multicohort)
        vcf_path = inputs.as_path(target=multicohort, stage=JointGenotyping, key='vcf')
        j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=get_batch().read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
                },
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(multicohort),
            output_summary_path=outputs['qc_summary'],
            output_detail_path=outputs['qc_detail'],
        )
        return self.make_outputs(multicohort, data=outputs, jobs=[j])


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
        if sequencing_group.participant_id not in config_retrieve(['validation', 'sample_map'], {}):
            return None

        return (
            get_multicohort().analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / 'hap.py'
            / f'{get_workflow().output_version}-{sequencing_group.id}.summary.csv'
        )

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        vcf_path = inputs.as_path(target=get_multicohort(), stage=JointGenotyping, key='vcf')

        jobs = happy(
            b=get_batch(),
            sequencing_group=sequencing_group,
            vcf_or_gvcf=get_batch().read_input_group(
                **{
                    'vcf.gz': str(vcf_path),
                    'vcf.gz.tbi': str(vcf_path) + '.tbi',
                },
            ),
            is_gvcf=False,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_path=self.expected_outputs(sequencing_group),
        )
        if not jobs:
            return self.make_outputs(sequencing_group)
        else:
            return self.make_outputs(sequencing_group, self.expected_outputs(sequencing_group), jobs)


@stage(required_stages=[JointVcfQC, JointVcfHappy], analysis_type='qc', analysis_keys=['json'])
class JointVcfMultiQC(MultiCohortStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        return {
            'html': self.web_prefix / 'multiqc.html',
            'json': self.prefix / 'multiqc_data.json',
            'checks': self.prefix / '.checks',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        outputs = self.expected_outputs(multicohort)
        if base_url := multicohort.analysis_dataset.web_url():
            html_url = str(outputs['html']).replace(str(self.web_prefix), base_url)
        else:
            html_url = None

        paths = []
        ending_to_trim = set()  # endings to trim to get sample names

        paths.append(inputs.as_path(multicohort, JointVcfQC, 'qc_detail'))

        for sequencing_group in multicohort.get_sequencing_groups():
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
            out_json_path=outputs['json'],
            out_html_path=outputs['html'],
            out_html_url=html_url,
            out_checks_path=outputs['checks'],
            job_attrs=self.get_job_attrs(multicohort),
            sequencing_group_id_map=multicohort.rich_id_map(),
            extra_config={'table_columns_visible': {'Picard': True}},
            dataset=multicohort.analysis_dataset,
            label='Joint VCF',
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)
