#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

import logging
from os.path import join
from typing import List

import click
from cloudpathlib import CloudPath

from cpg_pipes import buckets
from cpg_pipes.jobs import pedigree
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.stage import DatasetStage, StageInput, StageOutput
from cpg_pipes.pipeline.pipeline import stage, Pipeline
from cpg_pipes.pipeline.cli_opts import pipeline_click_options

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class CramPedCheckStage(DatasetStage):
    def expected_result(self, dataset: Dataset):
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        path_by_sid = dict()
        for s in dataset.get_samples():
            path = CloudPath(f'gs://cpg-{dataset.name}-main/cram/{s.id}.somalier')
            if buckets.file_exists(path):
                path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            dataset,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(CRAMs)',
            ignore_missing=self.pipe.skip_samples_with_missing_input,
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=somalier_samples_path, jobs=[j])


@stage
class GvcfPedCheckStage(DatasetStage):
    def expected_result(self, dataset: Dataset):
        pass

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        path_by_sid = dict()
        for s in dataset.get_samples():
            path = CloudPath(f'gs://cpg-{dataset.name}-main/gvcf/{s.id}.somalier')
            if buckets.file_exists(path):
                path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            dataset,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(GVCFs)',
            ignore_missing=self.pipe.skip_samples_with_missing_input,
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=somalier_samples_path, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    input_datasets: List[str],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_datasets

    title = f'Pedigree checks for {", ".join(input_datasets)}'

    Pipeline(
        name='pedigree_check',
        description=title,
        input_datasets=input_datasets,
        output_version=output_version,
        stages_in_order=[
            CramPedCheckStage,
            GvcfPedCheckStage,
        ],
        **kwargs,
    ).submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
