#!/usr/bin/env python3

"""
Batch pipeline check pedigree
"""

import logging
from os.path import join
from typing import Optional, List, Collection, Tuple, Dict, Union

import click
from hailtop.batch.job import Job

from cpg_pipes import utils, resources
from cpg_pipes.jobs import pedigree, wrap_command
from cpg_pipes.pipeline import Namespace, Pipeline, Project, \
    ProjectStage, CohortStage

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
        dep_paths_by_stage: Dict[str, Dict[str, str]] = None,
        dep_jobs: Optional[List[str]] = None,
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
        dep_paths_by_stage: Dict[str, Dict[str, str]] = None,
        dep_jobs: Optional[List[str]] = None,
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
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice([n.lower() for n in Namespace.__members__]),
    callback=lambda c, p, v: getattr(Namespace, v.upper()) if v else None,
    help='The bucket namespace to write the results to',
)
@click.option(
    '--analysis-project',
    'analysis_project',
    default='seqr',
    help='SM project name to write the intermediate/joint-calling analysis entries to',
)
@click.option(
    '--input-project',
    'input_projects',
    multiple=True,
    required=True,
    help='Only read samples that belong to the project(s). Can be set multiple times.',
)
@click.option(
    '--output-project',
    'output_projects',
    multiple=True,
    help='Only create ES indicies for the project(s). Can be set multiple times. '
    'Defaults to --input-projects. The name of the ES index will be suffixed '
    'with the dataset version (set by --version)',
)
@click.option(
    '--first-stage',
    'first_stage',
    help='Only pick results from the previous stages if they exist. '
    'If not, skip such samples',
)
@click.option(
    '--last-stage',
    'last_stage',
    help='Finish the pipeline after this stage',
)
@click.option(
    '--skip-sample',
    '-S',
    'skip_samples',
    multiple=True,
    help='Don\'t process specified samples. Can be set multiple times.',
)
@click.option(
    '--force-sample',
    'force_samples',
    multiple=True,
    help='Force reprocessing these samples. Can be set multiple times.',
)
@click.option(
    '--output-version',
    'output_version',
    type=str,
    default='v0',
    help='Suffix the outputs with this version tag. Useful for testing',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--ped-file',
    'ped_files',
    multiple=True,
)
@click.option(
    '--check-smdb-seq-existence/--no-check-smdb-seq-existence',
    'check_smdb_seq_existence',
    default=False,
    is_flag=True,
    help='Check that files in sequence.meta exist'
)
@click.option(
    '--skip-samples-without-seq-input',
    'skip_samples_without_seq_input',
    default=False,
    is_flag=True,
    help='If sequence.meta files for a sample don\'t exist, remove this sample '
         'instead of failing'
)
@click.option(
    '--check-intermediate-existence/--no-check-intermediate-existence',
    'check_intermediate_existence',
    default=True,
    is_flag=True,
    help='Before running a job, check for an intermediate output before submitting it, '
         'and if it exists on a bucket, submit a [reuse] job instead. Works well with '
         '--previous-batch-tsv/--previous-batch-id options.',
)
@click.option(
    '--update-smdb-analyses/--no-update-smdb-analyses',
    'update_smdb_analyses',
    is_flag=True,
    default=True,
    help='Create analysis entries for queued/running/completed jobs'
)
@click.option(
    '--validate-smdb-analyses/--no-validate-smdb-analyses',
    'validate_smdb_analyses',
    is_flag=True,
    default=False,
    help='Validate existing analysis entries by checking if a.output exists on '
         'the bucket. Set the analysis entry to "failure" if output doesn\'t exist'
)
def main(
    output_namespace: Namespace,
    analysis_project: str,
    input_projects: Collection[str],
    output_projects: Optional[Collection[str]],
    first_stage: Optional[str],
    last_stage: Optional[str],
    skip_samples: Collection[str],
    force_samples: Collection[str],
    output_version: str,
    keep_scratch: bool,
    dry_run: bool,
    ped_files: List[str],
    check_smdb_seq_existence: bool,
    skip_samples_without_seq_input: bool,
    check_intermediate_existence: bool,
    update_smdb_analyses: bool,
    validate_smdb_analyses: bool,
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
        analysis_project=analysis_project,
        name='pedigree_check_pipeline',
        output_version=output_version,
        namespace=output_namespace,
        keep_scratch=keep_scratch,
        title=title,
        smdb_update_analyses=update_smdb_analyses,
        smdb_check_seq_existence=check_smdb_seq_existence,
        skip_samples_without_seq_input=skip_samples_without_seq_input,
        validate_smdb_analyses=validate_smdb_analyses,
        check_intermediate_existence=check_intermediate_existence,
        first_stage=first_stage,
        last_stage=last_stage,
        input_projects=input_projects,
        skip_samples=skip_samples,
        force_samples=force_samples,
        ped_files=ped_files,
    )

    pipeline.run(dry_run)


class PedigreePipeline(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fingerprints_bucket = f'{self.analysis_bucket}/fingerprints'

        self.add_stages([
            CramPedCheckStage(self),
            GvcfPedCheckStage(self),
        ])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
