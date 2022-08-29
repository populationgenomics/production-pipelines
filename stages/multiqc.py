"""
Stage that summarises QC.
"""

import logging

from cpg_utils import Path
from cpg_utils.workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    StageInputNotFoundError,
    StageDecorator,
    DatasetStage,
    Dataset,
    exists,
)

from jobs.multiqc import multiqc
from .align import Align, qc_functions


@stage(
    required_stages=[
        Align,
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
        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        checks_path = self.expected_outputs(dataset)['checks']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        paths = []
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
                    if not exists(path):
                        logging.warning(
                            f'Output for stage {st.__name__} for {sample} '
                            f'does not exist: {path}'
                        )
                    else:
                        paths.append(path)
                        ending_to_trim.add(path.name.replace(sample.id, ''))

        assert ending_to_trim

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
