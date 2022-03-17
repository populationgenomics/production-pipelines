#!/usr/bin/env python3

"""
Batch pipeline to run WGS QC.
"""

import logging

import click
from cloudpathlib import CloudPath

from cpg_pipes.buckets import exists
from cpg_pipes.jobs.multiqc import multiqc
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.pipeline import stage, Pipeline, DatasetStage
from cpg_pipes.pipeline.stage import StageInput, StageOutput
from cpg_pipes.stages.cramqc import SamtoolsStats, PicardWgsMetrics, VerifyBamId
from cpg_pipes.stages.fastqc import FastQC
from pipelines.somalier import CramSomalierPedigree, CramSomalierAncestry

logger = logging.getLogger(__file__)


@stage(
    required_stages=[
        FastQC,
        SamtoolsStats, 
        PicardWgsMetrics,
        VerifyBamId,
        CramSomalierPedigree, 
        CramSomalierAncestry,
    ], 
    forced=True
)
class MultiQC(DatasetStage):
    """
    Run MultiQC to summarise all QC.
    """
    def expected_result(self, dataset: Dataset) -> dict[str, CloudPath]:
        """
        Expected to produce an HTML and a correponding JSON file.
        """
        return {
            'html': dataset.get_web_bucket() / 'qc' / 'multiqc.html',
            'json': dataset.get_analysis_bucket() / 'qc' / 'multiqc_data.json'
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Call a function from the `jobs` module using inputs from `cramqc` 
        and `somalier` stages.
        """
        somalier_samples = inputs.as_path(dataset, CramSomalierPedigree, id='samples')
        somalier_pairs = inputs.as_path(dataset, CramSomalierPedigree, id='pairs')
        somalier_ancestry = inputs.as_path(dataset, CramSomalierAncestry, id='tsv')

        json_path = self.expected_result(dataset)['json']
        html_path = self.expected_result(dataset)['html']
        if base_url := dataset.get_web_url():
            html_url = str(html_path).replace(str(dataset.get_web_bucket()), base_url)
        else:
            html_url = None

        paths = [somalier_samples, somalier_pairs, somalier_ancestry]
        ending_to_trim = set()  # endings to trim to get sample names
        for sample in dataset.get_samples():
            for path in [
                inputs.as_path(sample, FastQC, id='zip'),
                inputs.as_path(sample, SamtoolsStats),
                inputs.as_path(sample, PicardWgsMetrics),
                inputs.as_path(sample, VerifyBamId),
            ]:
                paths.append(path)
                ending_to_trim.add(path.name.replace(sample.id, ''))
        assert ending_to_trim
        modules_to_trim_endings = {
            'fastqc/zip', 'samtools', 'picard/wgs_metrics', 'verifybamid/selfsm'
        }

        j = multiqc(
            self.b,
            dataset_name=dataset.name,
            tmp_bucket=dataset.get_tmp_bucket(),
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
        )
        
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):
    """
    Entry point, decorated by pipeline click options.
    """
    pipeline = Pipeline(
        name='cram_qc',
        description='CRAM QC',
        **kwargs,
    )
    if pipeline.skip_samples_with_missing_input:
        for sample in pipeline.get_all_samples():
            if not exists(sample.get_cram_path().path):
                logger.warning(
                    f'Could not find CRAM, skipping sample {sample.id}'
                )
                sample.active = False

    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
