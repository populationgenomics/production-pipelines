#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples
"""

import logging
from os.path import join
from typing import Optional, Collection

import click

from cpg_pipes import utils
from cpg_pipes.jobs import pedigree
from cpg_pipes.pipeline import Project, \
    ProjectStage, pipeline_click_options, stage, StageResults, run_pipeline

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class CramPedCheckStage(ProjectStage):
    def get_expected_output(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageResults) -> StageResults:
        # path_by_sid = dict()
        # for s in project.samples:
        #     path = f'gs://cpg-{project.name}-main/cram/{s.id}.cram'
        #     if utils.file_exists(path):
        #         path_by_sid[s.id] = path
        path_by_sid = dict()
        for s in project.samples:
            path = f'gs://cpg-{project.name}-main/cram/{s.id}.somalier'
            # if utils.file_exists(path):
            path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(CRAMs)',
            ignore_missing=self.pipe.skip_samples_without_seq_input,
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


@stage
class GvcfPedCheckStage(ProjectStage):
    def get_expected_output(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageResults) -> StageResults:
        path_by_sid = dict()
        for s in project.samples:
            path = f'gs://cpg-{project.name}-main/gvcf/{s.id}.g.vcf.gz'
            if utils.file_exists(path):
                path_by_sid[s.id] = path

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(GVCFs)',
            ignore_missing=self.pipe.skip_samples_without_seq_input,
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    input_projects: Collection[str],
    output_projects: Optional[Collection[str]],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_projects

    title = f'Pedigree checks for {", ".join(input_projects)}'

    if not output_projects:
        output_projects = input_projects
    if not all(op in input_projects for op in output_projects):
        logger.critical(
            f'All output projects must be contained within the specified input '
            f'projects. Input project: {input_projects}, output projects: '
            f'{output_projects}'
        )

    run_pipeline(
        name='pedigree_check',
        title=title,
        input_projects=input_projects,
        output_version=output_version,
        stages_in_order=[
            CramPedCheckStage,
            GvcfPedCheckStage,
        ],
        **kwargs,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
