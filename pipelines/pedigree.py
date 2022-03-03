#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

import logging
from os.path import join
from typing import List

import click

from cpg_pipes import buckets
from cpg_pipes.jobs import pedigree
from cpg_pipes.pipeline.project import Project
from cpg_pipes.pipeline.stage import ProjectStage, StageInput, StageOutput
from cpg_pipes.pipeline.pipeline import stage, Pipeline
from cpg_pipes.pipeline.cli_opts import pipeline_click_options

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class CramPedCheckStage(ProjectStage):
    def expected_result(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        path_by_sid = dict()
        for s in project.get_samples():
            path = f'gs://cpg-{project.name}-main/cram/{s.id}.somalier'
            if buckets.file_exists(path):
                path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(CRAMs)',
            ignore_missing=self.pipe.skip_missing_input,
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


@stage
class GvcfPedCheckStage(ProjectStage):
    def expected_result(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        path_by_sid = dict()
        for s in project.get_samples():
            path = f'gs://cpg-{project.name}-main/gvcf/{s.id}.somalier'
            if buckets.file_exists(path):
                path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediates,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(GVCFs)',
            ignore_missing=self.pipe.skip_missing_input,
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    input_projects: List[str],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_projects

    title = f'Pedigree checks for {", ".join(input_projects)}'

    Pipeline(
        name='pedigree_check',
        description=title,
        input_projects=input_projects,
        output_version=output_version,
        stages_in_order=[
            CramPedCheckStage,
            GvcfPedCheckStage,
        ],
        **kwargs,
    ).submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
