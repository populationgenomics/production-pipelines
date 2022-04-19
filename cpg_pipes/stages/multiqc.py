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
from cpg_pipes.stages.cram_qc import SamtoolsStats, PicardWgsMetrics, VerifyBamId
from cpg_pipes.stages.fastqc import FastQC
from cpg_pipes.targets import Dataset
from pipelines.somalier import CramSomalierPedigree

logger = logging.getLogger(__file__)


@stage(
    required_stages=[
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
            'html': dataset.get_web_bucket() / 'qc' / 'multiqc.html',
            'json': dataset.get_bucket() / 'qc' / 'multiqc_data.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Call a function from the `jobs` module using inputs from `cramqc`
        and `somalier` stages.
        """
        somalier_samples = inputs.as_path(dataset, CramSomalierPedigree, id='samples')
        somalier_pairs = inputs.as_path(dataset, CramSomalierPedigree, id='pairs')

        json_path = self.expected_outputs(dataset)['json']
        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.get_web_url():
            html_url = str(html_path).replace(str(dataset.get_web_bucket()), base_url)
        else:
            html_url = None

        paths = [
            somalier_samples,
            somalier_pairs,
        ]
        ending_to_trim = set()  # endings to trim to get sample names
        for sample in dataset.get_samples():
            for st, key in [
                (SamtoolsStats, None),
                (PicardWgsMetrics, None),
                (VerifyBamId, None),
                (FastQC, 'zip'),
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
            'samtools',
            'picard/wgs_metrics',
            'verifybamid/selfsm',
        }
        
        # Building sample maps to MultiQC bulk rename. Only adding samples to the map
        # if the extrenal/participant IDs are differrent:
        external_id_map = {
            s.id: s.external_id for s in dataset.get_samples()
            if s.id != s.external_id
        }
        participant_id_map = {
            s.id: s.participant_id 
            for s in dataset.get_samples()
            if s.id != s.participant_id and s.participant_id != s.external_id
        }
        sid_maps = dict()
        if external_id_map:
            sid_maps['External ID'] = external_id_map
        if participant_id_map:
            sid_maps['Participant ID'] = participant_id_map

        j = multiqc(
            self.b,
            tmp_bucket=dataset.get_tmp_bucket(),
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            dataset_name=dataset.name,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(dataset),
            status_reporter=self.status_reporter,
            sample_id_maps=sid_maps,
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[j])
