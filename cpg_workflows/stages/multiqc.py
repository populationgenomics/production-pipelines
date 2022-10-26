"""
Stage that summarises QC.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    StageInputNotFoundError,
    StageDecorator,
    DatasetStage,
    Dataset,
    Cohort,
)

from cpg_workflows.jobs.multiqc import multiqc
from .align import Align, qc_functions
from .genotype import Genotype, GvcfHappy
from .joint_genotyping import JointGenotyping, JointVcfHappy
from .somalier import SomalierPedigree


@stage(
    required_stages=[
        Align,
        SomalierPedigree,
    ],
    forced=True,
)
class CramMultiQC(DatasetStage):
    """
    Run MultiQC to summarise all CRAM QC.
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
            stage_by_key: dict[str, StageDecorator] = {}
            for qc in qc_functions():
                for key, out in qc.outs.items():
                    if out:
                        stage_by_key[key] = Align
                        modules_to_trim_endings.add(out.multiqc_key)

            for key, st in stage_by_key.items():
                try:
                    path = inputs.as_path(sample, st, key)
                except StageInputNotFoundError:  # allow missing inputs
                    logging.warning(
                        f'Output for stage {st.__name__} not found for {sample}, '
                        f'skipping'
                    )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sample.id, ''))

        assert ending_to_trim

        jobs = multiqc(
            self.b,
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
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(
    required_stages=[
        Genotype,
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
            for st, key in [
                (Genotype, 'qc_detail'),
                (GvcfHappy, None),
            ]:
                try:
                    path = inputs.as_path(sample, st, key)
                except StageInputNotFoundError:  # allow missing inputs
                    if st != GvcfHappy:
                        logging.warning(
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

        jobs = multiqc(
            self.b,
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


@stage(
    required_stages=[
        JointGenotyping,
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
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}

        h = cohort.alignment_inputs_hash()
        return {
            'html': cohort.analysis_dataset.web_prefix() / 'qc' / 'jc' / 'multiqc.html',
            'json': cohort.analysis_dataset.prefix()
            / 'qc'
            / 'jc'
            / h
            / 'multiqc_data.json',
            'checks': cohort.analysis_dataset.prefix() / 'qc' / 'jc' / h / '.checks',
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

        paths.append(inputs.as_path(cohort, JointGenotyping, 'qc_detail'))

        for sample in cohort.get_samples():
            try:
                path = inputs.as_path(sample, JointVcfHappy)
            except StageInputNotFoundError:
                pass
            else:
                paths.append(path)
                ending_to_trim.add(path.name.replace(sample.id, ''))

        jobs = multiqc(
            self.b,
            tmp_prefix=self.tmp_prefix,
            paths=paths,
            ending_to_trim=ending_to_trim,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=checks_path,
            job_attrs=self.get_job_attrs(cohort),
            sample_id_map=cohort.rich_id_map(),
            extra_config={'table_columns_visible': {'Picard': True}},
            dataset=cohort.analysis_dataset,
            label='Joint VCF',
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
