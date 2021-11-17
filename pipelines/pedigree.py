#!/usr/bin/env python3

"""
Batch pipeline check pedigree
"""

import logging
from typing import Optional, List, Collection, Tuple, Dict, Union, Any

import click
from hailtop.batch.job import Job

from cpg_pipes import utils, resources
from cpg_pipes.jobs import pedigree, wrap_command
from cpg_pipes.pipeline import Namespace, Pipeline, Project, \
    ProjectStage, CohortStage, pipeline_click_options

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class CramPedCheckStage(ProjectStage):
    # defining __init__ only to make sure that .pipe instance is of type
    # PedigreePipeline and not Pipeline, so it has .fingerprints_bucket
    def __init__(self, pipe: 'PedigreePipeline'):
        super().__init__(pipe)
        self.pipe = pipe

    def get_expected_output(self, *args):
        pass

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[Any, Dict[str, str]] = None,
        dep_jobs: Optional[List[Job]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

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
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=dep_jobs or [],
            label='(CRAMs)',
            ignore_missing=self.pipe.skip_samples_without_seq_input,
        )
        return somalier_samples_path, [j]



class GvcfPedCheckStage(ProjectStage):
    # defining __init__ only to make sure that .pipe instance is of type
    # PedigreePipeline and not Pipeline, so it has .fingerprints_bucket
    def __init__(self, pipe: 'PedigreePipeline'):
        super().__init__(pipe)
        self.pipe = pipe

    def get_expected_output(self, *args):
        pass

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[Any, Dict[str, str]] = None,
        dep_jobs: Optional[List[Job]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

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
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=dep_jobs or [],
            label='(GVCFs)',
            ignore_missing=self.pipe.skip_samples_without_seq_input,
        )
        return somalier_samples_path, [j]


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

    pipeline = PedigreePipeline(
        name='seqr_loader',
        title=title,
        input_projects=input_projects,
        output_version=output_version,
        **kwargs,
    )
    pipeline.run()


class PedigreePipeline(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fingerprints_bucket = f'{self.analysis_bucket}/fingerprints'

        self.set_stages([
            CramPedCheckStage(self),
            GvcfPedCheckStage(self),
        ])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
