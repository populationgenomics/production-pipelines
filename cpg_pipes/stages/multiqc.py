"""
Stage that summarises QC.
"""
import logging

from cpg_utils.config import get_config

from cpg_pipes import Path
from cpg_pipes.jobs.multiqc import multiqc
from cpg_pipes.pipeline import (
    stage,
    StageInput,
    StageOutput,
    DatasetStage,
    CohortStage,
)
from cpg_pipes.pipeline.exceptions import StageInputNotFound
from cpg_pipes.stages.align import Align
from cpg_pipes.stages.cram_qc import SamtoolsStats, PicardWgsMetrics, VerifyBamId
from cpg_pipes.stages.fastqc import FastQC
from cpg_pipes.stages.somalier import (
    CramSomalierPedigree,
    GvcfSomalierPedigree,
)
from cpg_pipes.stages.variantqc import GvcfQc, JointVcfQc, GvcfHappy, JointVcfHappy
from cpg_pipes.targets import Dataset, Cohort

logger = logging.getLogger(__file__)


@stage(
    required_stages=[
        Align,
        FastQC,
        CramSomalierPedigree,
    ],
    forced=True,
)
class CramMultiQC(DatasetStage):
    """
    Run MultiQC to summarise all CRAM QC.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file. We write both to final
        "latest" location, and to a location parametrised by inputs hash, in order
        to enable reruns/reuses.
        """
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
        somalier_samples = inputs.as_path(dataset, CramSomalierPedigree, id='samples')
        somalier_pairs = inputs.as_path(dataset, CramSomalierPedigree, id='pairs')

        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        checks_path = self.expected_outputs(dataset)['checks']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        paths = [
            somalier_samples,
            somalier_pairs,
        ]
        ending_to_trim = set()  # endings to trim to get sample names
        for sample in dataset.get_samples():
            for st, key in [
                (FastQC, 'zip'),
                (Align, 'markduplicates_metrics'),
                (Align, 'samtools_stats'),
                (Align, 'picard_wgs_metrics'),
                (Align, 'verify_bamid'),
            ]:
                try:
                    path = inputs.as_path(sample, st, key)
                except StageInputNotFound:  # allow missing inputs
                    logger.warning(
                        f'Output for stage {st.__name__} not found for {sample}, '
                        f'skipping'
                    )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sample.id, ''))

        assert ending_to_trim
        modules_to_trim_endings = {
            'fastqc/zip',
            'picard/markdups',
            'samtools',
            'picard/wgs_metrics',
            'verifybamid/selfsm',
        }

        # Building sample map to MultiQC bulk rename. Only extending IDs if the
        # extrenal/participant IDs are differrent:
        jobs = multiqc(
            self.b,
            tmp_prefix=dataset.tmp_prefix(),
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
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(
    required_stages=[
        GvcfQc,
        GvcfSomalierPedigree,
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
        return {
            'html': dataset.web_prefix() / 'qc' / 'gvcf' / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'gvcf' / 'multiqc_data.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        somalier_samples = inputs.as_path(dataset, GvcfSomalierPedigree, id='samples')
        somalier_pairs = inputs.as_path(dataset, GvcfSomalierPedigree, id='pairs')

        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        paths = [
            somalier_samples,
            somalier_pairs,
        ]
        ending_to_trim = set()  # endings to trim to get sample names
        for sample in dataset.get_samples():
            for st, key in [
                (GvcfQc, 'detail'),
                (GvcfHappy, None),
            ]:
                try:
                    path = inputs.as_path(sample, st, key)
                except StageInputNotFound:  # allow missing inputs
                    if st != GvcfHappy:
                        logger.warning(
                            f'Output for stage {st.__name__} not found for {sample}, '
                            f'skipping'
                        )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sample.id, ''))

        assert ending_to_trim
        modules_to_trim_endings = {
            'picard/variant_calling_metrics',
            'happy',
        }

        # Building sample map to MultiQC bulk rename. Only extending IDs if the
        # extrenal/participant IDs are differrent:
        jobs = multiqc(
            self.b,
            tmp_prefix=dataset.tmp_prefix(),
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            dataset_name=dataset.name,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(dataset),
            sample_id_map=dataset.rich_id_map(),
            extra_config={'table_columns_visible': {'Picard': True}},
            name='GVCF QC report',
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(
    required_stages=[
        JointVcfQc,
        JointVcfHappy,
    ],
    forced=True,
)
class JointVcfMultiQC(CohortStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        h = cohort.alignment_inputs_hash()
        prefix = cohort.analysis_dataset.prefix() / 'qc' / 'jc'
        web_prefix = cohort.analysis_dataset.web_prefix() / 'qc' / 'jc'
        return {
            'html': web_prefix / 'multiqc.html',
            'json': prefix / 'multiqc_data.json',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        json_path = self.expected_outputs(cohort)['json']
        html_path = self.expected_outputs(cohort)['html']
        if base_url := cohort.analysis_dataset.web_url():
            html_url = str(html_path).replace(
                str(cohort.analysis_dataset.web_prefix()), base_url
            )
        else:
            html_url = None

        paths = [
            inputs.as_path(cohort, JointVcfQc, 'detail'),
            inputs.as_path(cohort, JointVcfHappy),
        ]

        # Building sample map to MultiQC bulk rename. Only extending IDs if the
        # extrenal/participant IDs are differrent:
        jobs = multiqc(
            self.b,
            tmp_prefix=self.tmp_prefix(),
            paths=paths,
            dataset_name=cohort.name,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(cohort),
            sample_id_map=cohort.rich_id_map(),
            name='Joint VCF QC report',
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
