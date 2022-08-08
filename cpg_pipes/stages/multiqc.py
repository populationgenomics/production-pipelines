"""
Stage that summarises QC.
"""

import logging

from cpg_pipes import Path
from cpg_pipes.jobs.multiqc import multiqc
from cpg_pipes.pipeline import (
    stage,
    StageInput,
    StageOutput,
    DatasetStage,
)
from cpg_pipes.pipeline.exceptions import StageInputNotFound
from cpg_pipes.stages.align import Align
from cpg_pipes.stages.cram_qc import SamtoolsStats, PicardWgsMetrics, VerifyBamId
from cpg_pipes.stages.fastqc import FastQC
from cpg_pipes.stages.somalier import CramSomalierPedigree
from cpg_pipes.targets import Dataset

logger = logging.getLogger(__file__)


@stage(
    required_stages=[
        Align,
        FastQC,
        SamtoolsStats,
        PicardWgsMetrics,
        VerifyBamId,
        CramSomalierPedigree,
    ],
    forced=True,
)
class MultiQC(DatasetStage):
    """
    Run MultiQC to summarise all QC.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        return {
            'html': dataset.web_prefix() / 'qc' / 'multiqc.html',
            'json': dataset.prefix() / 'qc' / 'multiqc_data.json',
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
                (SamtoolsStats, None),
                (PicardWgsMetrics, None),
                (VerifyBamId, None),
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
        j = multiqc(
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
            status_reporter=self.status_reporter,
            sample_id_map=dataset.rich_id_map(),
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[j])
