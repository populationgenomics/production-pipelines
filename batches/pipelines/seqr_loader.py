#!/usr/bin/env python3

"""
Batch pipeline to laod data info seqr
"""

import logging
import os
import sys
import time
from dataclasses import dataclass
from os.path import join
from typing import Optional, List, Collection, Tuple, Dict

import click
import hail as hl
from analysis_runner import dataproc
from hailtop.batch.job import Job

from cpg_pipes import utils, resources
from cpg_pipes.jobs import align, split_intervals, haplotype_caller, \
    pedigree
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.pipeline import Namespace, Pipeline, Sample, \
    SampleStage, Project, CohortStage, AnalysisType, ProjectStage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class CramStage(SampleStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            analysis_type=AnalysisType.CRAM,
        )

    def get_expected_output(self, sample: Sample):
        return f'{sample.project.get_bucket()}/cram/{sample.id}.cram'

    @staticmethod
    def get_fingerprint_output(output_path):
        return output_path.replace('.cram', '.somalier')

    def add_jobs(
        self, 
        sample: Sample,
        dep_paths_by_stage=None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        if not sample.seq_info:
            logger.critical(f'No sequence record for {sample.id}')
            sys.exit(1)
        
        alignment_input = sample.seq_info.parse_reads_from_metadata()
        if not alignment_input:
            if not self.pipe.skip_samples_without_seq_input:
                logger.critical(f'Could not find read data for sample {sample.id}')
                sys.exit(1)
            else:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.project.samples = [
                    s for s in sample.project.samples if s is not sample
                ]
                return None, []

        expected_path = self.get_expected_output(sample)
        cram_job = align.bwa(
            b=self.pipe.b,
            alignment_input=alignment_input,
            output_path=expected_path,
            sample_name=sample.id,
            project_name=sample.project.name,
            overwrite=not self.pipe.check_intermediate_existence,
            check_existence=False,
            smdb=self.pipe.db,
            prev_batch_jobs=self.pipe.prev_batch_jobs,
        )
        
        fingerprint_job, fingerprint_path = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=expected_path,
            overwrite=not self.pipe.check_intermediate_existence,
            label='(CRAMs)',
            depends_on=[cram_job],
        )
        
        return expected_path, [cram_job, fingerprint_job]


class CramPedCheckStage(ProjectStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[CramStage],
        )
        self.pipe = pipe

    def get_expected_output(self, project: Project):
        pass

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[str, Dict[str, str]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        cram_by_sid = dep_paths_by_stage[CramStage.get_name()]
        fingerprint_path_by_sid = {
            sid: CramStage.get_fingerprint_output(cram_path)
            for sid, cram_path in cram_by_sid
        }

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=fingerprint_path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=dep_jobs or [],
            label='(CRAMs)',
        )
        return somalier_samples_path, [j]


class GvcfStage(SampleStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            analysis_type=AnalysisType.GVCF,
            requires_stages=[CramStage],
        )

    def get_expected_output(self, sample: Optional[Sample]):
        return f'{sample.project.get_bucket()}/gvcf/{sample.id}.g.vcf.gz'

    @staticmethod
    def get_fingerprint_output(output_path):
        return output_path.replace('.g.vcf.gz', '.somalier')
    
    def add_jobs(
        self, 
        sample: Sample,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        cram_path = dep_paths_by_stage[CramStage.get_name()]

        if self.pipe.hc_intervals is None and self.pipe.config.hc_shards_num > 1:
            self.pipe.hc_intervals = split_intervals.get_intervals(
                b=self.pipe.b,
                scatter_count=self.pipe.config.hc_shards_num,
            )
        expected_path = self.get_expected_output(sample)
        gvcf_job = haplotype_caller.produce_gvcf(
            b=self.pipe.b,
            output_path=expected_path,
            sample_name=sample.id,
            project_name=sample.project.name,
            cram_path=cram_path,
            crai_path=cram_path + '.crai',
            intervals=self.pipe.hc_intervals,
            number_of_intervals=self.pipe.config.hc_shards_num,
            tmp_bucket=self.pipe.tmp_bucket,
            overwrite=not self.pipe.check_intermediate_existence,
            check_existence=False,
            depends_on=dep_jobs,
            smdb=self.pipe.db,
        )
        fingerprint_job, fingerprint_path = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=expected_path,
            overwrite=not self.pipe.check_intermediate_existence,
            label='(GVCFs)',
            depends_on=[gvcf_job],
        )
        return expected_path, [gvcf_job, fingerprint_job]


class GvcfPedCheckStage(ProjectStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[GvcfStage],
        )
        self.pipe = pipe

    def get_expected_output(self, *args):
        pass

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[str, Dict[str, str]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        gvcf_by_sid = dep_paths_by_stage[GvcfStage.get_name()]
        fingerprint_path_by_sid = {
            sid: CramStage.get_fingerprint_output(gvcf_path)
            for sid, gvcf_path in gvcf_by_sid
        }

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=fingerprint_path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=dep_jobs or [],
            label='(GVCFs)'
        )
        return somalier_samples_path, [j]


class JointGenotypingStage(CohortStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[GvcfStage],
            analysis_type=AnalysisType.JOINT_CALLING,
        )

    def get_expected_output(self, *args):
        samples_hash = utils.hash_sample_ids(self.pipe.get_all_sample_ids())
        expected_jc_vcf_path = f'{self.pipe.tmp_bucket}/joint_calling/{samples_hash}.vcf.gz'
        return expected_jc_vcf_path
    
    @staticmethod
    def make_expected_siteonly_output(output_path):
        return output_path.replace('.vcf.gz', '-siteonly.vcf.gz')

    def add_jobs(
        self,
        pipe: Pipeline,
        dep_paths_by_stage: Dict[str, Dict[str, str]] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        
        gvcf_path_by_sid = dep_paths_by_stage[GvcfStage.get_name()]

        not_found_gvcfs = []
        for sid, gvcf_path in gvcf_path_by_sid.items():
            if gvcf_path is None:
                logger.error(f'Joint genotyping: could not find GVCF for {sid}')
                not_found_gvcfs.append(sid)
        if not_found_gvcfs:
            logger.critical(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting')
            sys.exit(1)

        expected_path = self.get_expected_output(pipe)
        jc_job = make_joint_genotyping_jobs(
            b=self.pipe.b,
            out_vcf_path=expected_path,
            out_siteonly_vcf_path=self.make_expected_siteonly_output(expected_path),
            samples=self.pipe.get_all_samples(),
            genomicsdb_bucket=f'{self.pipe.analysis_bucket}/genomicsdbs',
            tmp_bucket=self.pipe.tmp_bucket,
            gvcf_by_sid=gvcf_path_by_sid,
            local_tmp_dir=self.pipe.local_tmp_dir,
            overwrite=not self.pipe.check_intermediate_existence,
            depends_on=dep_jobs or [],
            smdb=self.pipe.db,
            tool=JointGenotyperTool.GnarlyGenotyper 
            if self.pipe.config.use_gnarly 
            else JointGenotyperTool.GenotypeGVCFs,
        )
        return expected_path, [jc_job]


class VqsrStage(CohortStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[JointGenotypingStage],
        )

    def get_expected_output(self, *args):
        samples_hash = utils.hash_sample_ids(self.pipe.get_all_sample_ids())
        expected_jc_vcf_path = f'{self.pipe.tmp_bucket}/vqsr/{samples_hash}-site-only.vcf.gz'
        return expected_jc_vcf_path
    
    def add_jobs(
        self,
        pipe: Pipeline,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        jc_vcf_path = dep_paths_by_stage.get(JointGenotypingStage.get_name())
        siteonly_vcf_path = JointGenotypingStage.make_expected_siteonly_output(jc_vcf_path)

        tmp_vqsr_bucket = f'{self.pipe.tmp_bucket}/vqsr'
        logger.info(f'Queueing VQSR job')
        expected_path = self.get_expected_output(pipe)
        vqsr_job = make_vqsr_jobs(
            b=self.pipe.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            work_bucket=tmp_vqsr_bucket,
            web_bucket=tmp_vqsr_bucket,
            gvcf_count=len(self.pipe.get_all_samples()),
            depends_on=dep_jobs or [],
            scatter_count=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
            output_vcf_path=expected_path,
            use_as_annotations=self.pipe.config.use_as_vqsr,
            overwrite=not self.pipe.check_intermediate_existence,
        )
        return expected_path, [vqsr_job]


class AnnotateCohortStage(CohortStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[JointGenotypingStage, VqsrStage],
        )
        self.anno_tmp_bucket = f'{self.pipe.tmp_bucket}/mt'

    def get_expected_output(self, *args):
        return f'{self.anno_tmp_bucket}/combined.mt'

    def add_jobs(
        self,
        pipe: Pipeline,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        
        checkpoints_bucket = f'{self.anno_tmp_bucket}/checkpoints'
        
        vcf_path = dep_paths_by_stage.get(JointGenotypingStage.get_name())
        vqsr_vcf_path = dep_paths_by_stage.get(VqsrStage.get_name())

        expected_path = self.get_expected_output(pipe)
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "vcf_to_mt.py")} '
            f'--vcf-path {vcf_path} '
            f'--site-only-vqsr-vcf-path {vqsr_vcf_path} '
            f'--dest-mt-path {expected_path} '
            f'--bucket {checkpoints_bucket} '
            f'--disable-validation '
            f'--make-checkpoints '
            + ('--overwrite ' if not self.pipe.check_intermediate_existence else ''),
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=50,
            job_name='Make MT and annotate cohort',
            vep='GRCh38',
            depends_on=dep_jobs or [],
        )
        return expected_path, [j]


class AnnotateProjectStage(ProjectStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
            requires_stages=[AnnotateCohortStage],
        )

    def get_expected_output(self, project: Project):
        return f'{project.get_bucket()}/mt/{project.name}.mt'

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        
        if project.name not in self.pipe.config.output_projects:
            logger.info(
                f'Skipping annotating project {project.name} because it is not'
                f'in the --output-projects: {self.pipe.config.output_projects}'
            )
            return None, []
        
        annotated_mt_path = dep_paths_by_stage.get(AnnotateCohortStage.get_name())

        # Make a list of project samples to subset from the entire matrix table
        sample_ids = [s.id for s in project.samples]
        proj_tmp_bucket = project.get_tmp_bucket()
        subset_path = f'{proj_tmp_bucket}/seqr-samples.txt'
        with hl.hadoop_open(subset_path, 'w') as f:
            f.write('\n'.join(sample_ids))

        expected_path = self.get_expected_output(project)    
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_projectmt.py")} '
            f'--mt-path {annotated_mt_path} '
            f'--out-mt-path {expected_path} '
            f'--subset-tsv {subset_path}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=20,
            job_name=f'{project.name}: annotate project',
            depends_on=dep_jobs or [],
        )
        return expected_path, [j]


class LoadToEsStage(ProjectStage):
    def __init__(self, pipe: 'SeqrLoaderPipeline'):
        super().__init__(
            pipe, 
        )

    def get_expected_output(self, project: Project):
        return None

    def add_jobs(
        self,
        project: Project,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs: Optional[List[str]] = None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:

        if project.name not in self.pipe.config.output_projects:
            logger.info(
                f'Skipping loading project {project.name} because it is not'
                f'in the --output-projects: {self.pipe.config.output_projects}'
            )
            return None, []

        project_mt_path = project.output_by_stage.get(AnnotateProjectStage.get_name())
        dep_jobs = project.jobs_by_stage.get(AnnotateProjectStage.get_name())

        timestamp = time.strftime("%Y%m%d-%H%M%S")
        j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "projectmt_to_es.py")} '
            f'--mt-path {project_mt_path} '
            f'--es-index {project.name}-{self.pipe.output_version}-{timestamp} '
            f'--es-index-min-num-shards 1 '
            f'{"--prod" if self.pipe.namespace == Namespace.MAIN else ""}',
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=2,
            job_name=f'{project.name}: create ES index',
            depends_on=dep_jobs or [],
            scopes=['cloud-platform'],
        )
        
        return None, [j]


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
@click.option(
    '--keep-scratch/--remove-scratch', 
    'keep_scratch', 
    default=True,
    is_flag=True,
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--ped-file',
    'ped_files',
    multiple=True,
)
@click.option(
    '--skip-ped-checks',
    'skip_ped_checks',
    is_flag=True,
    help='Skip checking provided sex and pedigree against the inferred one',
)
@click.option(
    '--hc-shards-num',
    'hc_shards_num',
    type=click.INT,
    default=resources.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
    help='Number of intervals to devide the genome for gatk HaplotypeCaller',
)
@click.option(
    '--use-gnarly/--no-use-gnarly',
    'use_gnarly',
    default=False,
    is_flag=True,
    help='Use GnarlyGenotyper instead of GenotypeGVCFs',
)
@click.option(
    '--use-as-vqsr/--no-use-as-vqsr',
    'use_as_vqsr',
    default=True,
    is_flag=True,
    help='Use allele-specific annotations for VQSR',
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
@click.option(
    '--previous-batch-tsv',
    'previous_batch_tsv_path',
    help='A list of previous successful attempts from another batch, dumped from '
         'from the Batch database (the "jobs" table joined on "job_attributes"). '
         'If the intermediate output for a job exists in a previous attempt, '
         'it will be passed forward, and a [reuse] job will be submitted.'
)
@click.option(
    '--previous-batch-id',
    'previous_batch_id',
    help='6-letter ID of the previous successful batch (corresponds to the directory '
         'name in the batch logs. e.g. feb0e9 in '
         'gs://cpg-seqr-main-tmp/hail/batch/feb0e9'
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
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    check_smdb_seq_existence: bool,
    skip_samples_without_seq_input: bool,
    check_intermediate_existence: bool,
    update_smdb_analyses: bool,
    validate_smdb_analyses: bool,
    previous_batch_tsv_path: Optional[str],
    previous_batch_id: Optional[str],
):  # pylint: disable=missing-function-docstring
    # Determine bucket paths

    assert input_projects

    title = f'Seqr loading: joint call from: {", ".join(input_projects)}'
    if output_projects:
        title += f', ES index for: {", ".join(output_projects)}'
    title += f', version {output_version}'
    
    if not output_projects:
        output_projects = input_projects
    if not all(op in input_projects for op in output_projects):
        logger.critical(
            f'All output projects must be contained within the specified input '
            f'projects. Input project: {input_projects}, output projects: '
            f'{output_projects}'
        )

    pipeline = SeqrLoaderPipeline(
        analysis_project=analysis_project,
        name='seqr_loader',
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
        config=SeqrLoaderConfig(
            hc_shards_num=hc_shards_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
            skip_ped_checks=skip_ped_checks,
            output_projects=output_projects,
        ),
        input_projects=input_projects,
        skip_samples=skip_samples,
        force_samples=force_samples,
        ped_files=ped_files,
        previous_batch_tsv_path=previous_batch_tsv_path,
        previous_batch_id=previous_batch_id,
    )

    pipeline.run(dry_run)


@dataclass
class SeqrLoaderConfig:
    hc_shards_num: int
    use_gnarly: bool
    use_as_vqsr: bool
    skip_ped_checks: bool
    output_projects: List[str]


class SeqrLoaderPipeline(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fingerprints_bucket = f'{self.analysis_bucket}/fingerprints'
        
        self.cram_path_by_sid = dict()
        self.cram_job_by_sid = dict()
        self.gvcf_path_by_sid = dict()
        self.gvcf_job_by_sid = dict()
        self.joint_called_vcf_path = None
        self.vqsr_site_only_vcf_path = None
        self.jc_job = None
        self.annotated_mt_path = None
        self.annotate_job = None
        self.hc_intervals = None

        self.add_stages([
            CramStage(self),
            CramPedCheckStage(self),
            GvcfStage(self),
            GvcfPedCheckStage(self),
            JointGenotypingStage(self),
            VqsrStage(self),
            AnnotateCohortStage(self),
            AnnotateProjectStage(self),
            LoadToEsStage(self),
        ])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
