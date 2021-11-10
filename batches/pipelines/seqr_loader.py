#!/usr/bin/env python3

"""
Driver for loading data into SEQR for the CPG. See the README for more information.

- 2021/04/16 Michael Franklin and Vlad Savelyev
"""

import logging
import sys
import time
from dataclasses import dataclass
from enum import Enum
from os.path import join
from typing import Optional, List, Collection

import pandas as pd
import click
import hail as hl
from analysis_runner import dataproc
from hailtop.batch.job import Job

from cpg_production_pipelines import utils, resources
from cpg_production_pipelines.jobs import align, split_intervals, haplotype_caller, \
    pedigree
from cpg_production_pipelines.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_production_pipelines.pipeline import Namespace, Pipeline, Sample, \
    SampleStage, Project, CohortStage, AnalysisType, ProjectStage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Stage(Enum):
    INPUT = 1
    CRAM = 2
    GVCF = 3
    GVCF_PED_CHECK = 4
    JOINT_CALLING = 5
    ANNOTATE = 6
    LOAD_TO_ES = 7


@click.command()
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice([n.lower() for n in Namespace.__members__]),
    callback=lambda c, p, v: getattr(Namespace, v) if v else None,
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
    '--start-from-stage',
    'first_stage',
    type=click.Choice([n.lower() for n in Stage.__members__]),
    callback=lambda c, p, v: getattr(Stage, v) if v else None,
    help='Only pick results from the previous stages if they exist. '
    'If not, skip such samples',
)
@click.option(
    '--end-with-stage',
    'last_stage',
    type=click.Choice([n.lower() for n in Stage.__members__]),
    callback=lambda c, p, v: getattr(Stage, v) if v else None,
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
    '--output-version',
    'output_version',
    type=str,
    default='v0',
    help='Suffix the outputs with this version tag. Useful for testing',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate Hail data',
)
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
@click.option('--vep-block-size', 'vep_block_size', type=click.INT)
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
    '--check-inputs-existence/--skip-check-inputs-existence',
    'check_inputs_existence',
    default=True,
    is_flag=True,
)
@click.option(
    '--update-smdb/--skip-update-smdb',
    'update_smdb',
    default=True,
    is_flag=True,
    help='Create analysis entries for queued/running/completed jobs'
)
@click.option(
    '--validate-smdb-analyses',
    'validate_smdb_analyses',
    is_flag=True,
    help='Validate existing analysis entries. Set to failure if output doesn\'t exist'
)
def main(
    output_namespace: Namespace,
    analysis_project: str,
    input_projects: Collection[str],
    output_projects: Optional[Collection[str]],
    first_stage: Optional[Stage],
    last_stage: Optional[Stage],
    skip_samples: Collection[str],
    output_version: str,
    keep_scratch: bool,
    overwrite: bool,
    dry_run: bool,
    ped_files: List[str],
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    check_inputs_existence: bool,
    update_smdb: bool,
    validate_smdb_analyses: bool,
):  # pylint: disable=missing-function-docstring
    # Determine bucket paths

    assert input_projects
    if output_projects:
        if not all(op in input_projects for op in output_projects):
            logger.critical(
                'All output projects must be contained within '
                'the specified input projects'
            )

    if output_namespace != Namespace.MAIN:
        analysis_project = f'{analysis_project}-test'
        input_projects = [f'{p}-test' for p in input_projects]
        output_projects = [f'{p}-test' for p in output_projects]

    pipeline = SeqrLoaderPipeline(
        analysis_project=analysis_project,
        name='seqr_loader',
        output_version=output_version,
        namespace=output_namespace,
        keep_scratch=keep_scratch,
        title=(
            f'Seqr loading. '
            f'{", ".join(input_projects)} -> '
            f'{", ".join(output_projects)}, '
            f'v{output_version}'
        ),
        do_update_analyses=update_smdb,
        do_check_existence=check_inputs_existence,
        validate_smdb_analyses=validate_smdb_analyses,
        overwrite=overwrite,
        first_stage=(first_stage or list(Stage.__members__.keys())[0]).name,
        last_stage=(last_stage or list(Stage.__members__.keys())[-1]).name,
        config=SeqrLoaderConfig(
            input_projects=input_projects,
            skip_samples=skip_samples,
            hc_shards_num=hc_shards_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
            skip_ped_checks=skip_ped_checks,
        ),
    )
    
    pipeline.populate_samples(
        input_projects=input_projects,
        skip_samples=skip_samples,
        namespace=output_namespace,
        ped_files=ped_files,
    )

    pipeline.run(dry_run)


@dataclass
class SeqrLoaderConfig:
    input_projects: List[str]
    skip_samples: List[str]
    hc_shards_num: int
    use_gnarly: bool
    use_as_vqsr: bool
    skip_ped_checks: bool


class JointCalledOutputType(Enum):
    pre_vqsr_vcf = 1
    vqsr_site_only_vcf = 2
    cohort_mt = 3
    project_mt = 4


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

        self.add_stages({
            Stage.CRAM: CramStage(self),
            Stage.CRAM_PED_CHECK: CramPedCheckStage(self),
            Stage.GVCF: GvcfStage(self),
            Stage.GVCF_PED_CHECK: GvcfPedCheckStage(self),
            Stage.JOINT_CALLING: JointGenotypingStage(self),
            Stage.ANNOTATE: AnnotateStage(self),
            Stage.LOAD_TO_ES: LoadToEsStage(self),
        })
        
    # def _vqsr(self):
    #     output = self._find_joint_called_output()
    # 
    #     if not self.joint_called_vcf_path:
    #         logger.critical('Joint-caled VCF is not found')
    #     
    #     tmp_vqsr_bucket = f'{self.tmp_bucket}/vqsr'
    #     logger.info(f'Queueing VQSR job')
    #     vqsr_job = make_vqsr_jobs(
    #         b=self.b,
    #         input_vcf_or_mt_path=jc_vcf_path,
    #         work_bucket=tmp_vqsr_bucket,
    #         web_bucket=tmp_vqsr_bucket,
    #         intervals=intervals_j.intervals,
    #         gvcf_count=len(samples),
    #         depends_on=[final_gathered_vcf_job],
    #         scatter_count=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
    #         output_vcf_path=output_path,
    #         use_as_annotations=use_as_vqsr,
    #         overwrite=overwrite,
    #     )
    #     if smdb:
    #         last_j = smdb.add_running_and_completed_update_jobs(
    #             b=b,
    #             analysis_type='joint-calling',
    #             output_path=out_jc_vcf_path,
    #             sample_names=sample_ids,
    #             first_j=intervals_j,
    #             last_j=last_j,
    #             depends_on=depends_on,
    #         )


class CramStage(SampleStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.CRAM.name, analysis_type=AnalysisType.CRAM)
    
    def get_expected_output(self, sample: Optional[Sample] = None):
        assert sample
        return f'{sample.project.get_bucket()}/cram/{sample.id}.cram'
    
    def add_jobs(self, sample: Sample):
        found_path = self.find_stage_output(sample)
        if not found_path:
            if not sample.seq_info:
                logger.critical(f'No sequence record for {sample.id}')
                return
            alignment_input = sample.seq_info.parse_reads_from_metadata()
            if not alignment_input:
                logger.critical(f'Could not find read data for sample {sample.id}')
                return
            expected_path = self.get_expected_output(sample)
            cram_job = align.bwa(
                b=self.pipe.b,
                alignment_input=alignment_input,
                output_path=expected_path,
                sample_name=sample.id,
                project_name=sample.project.name,
                overwrite=self.pipe.overwrite,
            )
            sample.jobs_by_stage[self.name] = [cram_job]
            found_path = expected_path
        sample.output_by_stage[self.name] = found_path


class GvcfStage(SampleStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.GVCF.name, analysis_type=AnalysisType.GVCF)

    def get_expected_output(self, sample: Optional[Sample] = None):
        assert sample
        return f'{sample.project.get_bucket()}/gvcf/{sample.id}.g.vcf.gz'

    def add_jobs(self, sample: Sample):
        found_path = self.find_stage_output(sample)
        if not found_path:
            cram_path = sample.output_by_stage.get(Stage.CRAM.name)
            if not cram_path:
                logger.critical(
                    f'No CRAM found for {sample.project.name}/{sample.id}. '
                    f'Required for GVCF analysis'
                )
                sys.exit(1)

            if self.pipe.hc_intervals is None and self.pipe.config.hc_shards_num > 1:
                self.pipe.hc_intervals = split_intervals.get_intervals(
                    b=self.pipe.b,
                    scatter_count=self.pipe.config.hc_shards_num,
                )
            expected_path = self.get_expected_output(sample)
            cram_jobs = sample.jobs_by_stage.get(Stage.CRAM.name, [])
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
                overwrite=self.pipe.overwrite,
                depends_on=cram_jobs,
                smdb=self.pipe.db,
            )
            sample.jobs_by_stage[self.name] = [gvcf_job]
            found_path = expected_path
        sample.output_by_stage[self.name] = found_path


class CramPedCheckStage(ProjectStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.CRAM_PED_CHECK.name)
        self.pipe = pipe

    def get_expected_output(self):
        pass

    def add_jobs(self, project: Project):
        cram_jobs = [s.jobs_by_stage[Stage.CRAM.name] for s in project.samples]
        file_by_sid = {s.id: s.output_by_stage[Stage.CRAM.name] for s in project.samples}
        pedigree.job(
            self.pipe.b,
            project,
            file_by_sid=file_by_sid,
            overwrite=self.pipe.overwrite,
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=cram_jobs,
        )


class GvcfPedCheckStage(ProjectStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.GVCF_PED_CHECK.name)
        self.pipe = pipe

    def get_expected_output(self):
        pass

    def add_jobs(self, project: Project):
        cram_jobs = [s.jobs_by_stage[Stage.GVCF.name] for s in project.samples]
        file_by_sid = {s.id: s.output_by_stage[Stage.GVCF.name] for s in project.samples}
        pedigree.job(
            self.pipe.b,
            project,
            file_by_sid=file_by_sid,
            overwrite=self.pipe.overwrite,
            fingerprints_bucket=self.pipe.fingerprints_bucket,
            web_bucket=self.pipe.web_bucket,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=cram_jobs,
        )


class JointGenotypingStage(CohortStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.JOINT_CALLING.name, analysis_type=AnalysisType.JOINT_CALLING)

    def get_expected_output(self):
        samples_hash = utils.hash_sample_ids(self.pipe.get_all_sample_ids())
        expected_jc_vcf_path = f'{self.pipe.tmp_bucket}/joint_calling/{samples_hash}.vcf.gz'
        return expected_jc_vcf_path
    
    @staticmethod
    def make_expected_vqsr_site_only_vcf_path(vcf_path):
        # Analysis in SMDB doesn't support multiple 'output', so we trust 
        # that the VQSR site-only VCF exists if the standard VCF exists
        return vcf_path.replace('.vcf.gz', '-vqsr-siteonly.vcf.gz')

    def add_jobs(self):
        found_joint_vcf_path = self.find_stage_output(self.pipe)
        samples = self.pipe.get_all_samples()
        if not found_joint_vcf_path:
            gvcf_path_by_sid = {
                s.id: s.output_by_stage.get(Stage.GVCF.name) for s in samples
            }
            
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
    
            expected_path = self.get_expected_output()
            expected_vqsr_path = self.make_expected_vqsr_site_only_vcf_path(expected_path)
            gvcf_jobs = [s.jobs_by_stage.get(Stage.GVCF.name) for s in samples]
            gvcf_jobs = [j[0] for j in gvcf_jobs if j]

            jc_job = make_joint_genotyping_jobs(
                b=self.pipe.b,
                out_jc_vcf_path=expected_path,
                out_vqsr_site_only_vcf_path=expected_vqsr_path,
                samples=self.pipe.get_all_samples(),
                genomicsdb_bucket=f'{self.pipe.analysis_bucket}/genomicsdbs',
                tmp_bucket=self.pipe.tmp_bucket,
                gvcf_by_sid=gvcf_path_by_sid,
                local_tmp_dir=self.pipe.local_tmp_dir,
                overwrite=self.pipe.overwrite,
                depends_on=gvcf_jobs,
                smdb=self.pipe.db,
                tool=JointGenotyperTool.GnarlyGenotyper 
                if self.pipe.config.use_gnarly 
                else JointGenotyperTool.GenotypeGVCFs,
            )
            self.pipe.jobs_by_stage[self.name] = [jc_job]
            found_joint_vcf_path = expected_path

        self.pipe.output_by_stage[self.name] = found_joint_vcf_path


class AnnotateStage(CohortStage):
    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, name=Stage.ANNOTATE.name)
        self.anno_tmp_bucket = f'{self.pipe.tmp_bucket}/mt'

    def get_expected_output(self, *args):
        return f'{self.anno_tmp_bucket}/combined.mt'

    def add_jobs(self):
        annotated_mt_path = self.find_stage_output(self.pipe)
        job_name = 'Make MT and annotate [reuse]'
        if self.pipe.can_reuse(annotated_mt_path):
            annotate_combined_j = self.pipe.b.new_job(job_name)
        else:
            checkpoints_bucket = f'{self.anno_tmp_bucket}/checkpoints'
            vcf_path = self.pipe.output_by_stage.get(Stage.JOINT_CALLING.name)
            if not self.pipe.can_reuse(vcf_path):
                logger.critical(
                    f'Could not find joint-called VCF path, cannot run annotation')
                sys.exit(1)
    
            vqsr_vcf_path = JointGenotypingStage.make_expected_vqsr_site_only_vcf_path(vcf_path)
            if not self.pipe.can_reuse(vqsr_vcf_path):
                logger.critical(
                    f'Could not find joint-called VCF path, cannot run annotation')
                sys.exit(1)
    
            annotate_combined_j = dataproc.hail_dataproc_job(
                self.pipe.b,
                f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "vcf_to_mt.py")} '
                f'--vcf-path {vcf_path} '
                f'--site-only-vqsr-vcf-path {vqsr_vcf_path} '
                f'--dest-mt-path {annotated_mt_path} '
                f'--bucket {checkpoints_bucket} '
                '--disable-validation '
                '--make-checkpoints '
                + ('--overwrite ' if self.pipe.overwrite else ''),
                max_age='16h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=50,
                job_name=job_name,
                vep='GRCh38',
                depends_on=self.pipe.jobs_by_stage.get(Stage.JOINT_CALLING.name, []),
            )
        self.pipe.output_by_stage[self.name] = annotated_mt_path
        self.pipe.jobs_by_stage[self.name] = [annotate_combined_j]


class LoadToEsStage(ProjectStage):
    def get_expected_output(self, *args):
        return None

    def __init__(self, pipe: SeqrLoaderPipeline):
        super().__init__(pipe, Stage.LOAD_TO_ES.name)

    def add_jobs(self, project: Project):
        annotated_mt_path = self.pipe.output_by_stage.get(Stage.ANNOTATE.name)
        proj_name = project.name
        proj_name = proj_name if not project.is_test else f'{proj_name}-test'
        proj_bucket = project.get_bucket()
        proj_tmp_bucket = project.get_tmp_bucket()
        project_mt_path = f'{proj_bucket}/mt/{proj_name}.mt'
        # Make a list of project samples to subset from the entire matrix table
        sample_ids = [s.id for s in project.samples]
        subset_path = f'{proj_tmp_bucket}/seqr-samples.txt'
        with hl.hadoop_open(subset_path, 'w') as f:
            f.write('\n'.join(sample_ids))

        annotate_project_j = dataproc.hail_dataproc_job(
            self.pipe.b,
            f'batch_seqr_loader/scripts/mt_to_projectmt.py '
            f'--mt-path {annotated_mt_path} '
            f'--out-mt-path {project_mt_path}'
            f'--subset-tsv {subset_path}',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            job_name=f'{proj_name}: annotate project',
            depends_on=self.pipe.jobs_by_stage.get(Stage.ANNOTATE, []),
        )

        timestamp = time.strftime("%Y%m%d-%H%M%S")
        dataproc.hail_dataproc_job(
            self.pipe.b,
            f'batch_seqr_loader/scripts/projectmt_to_es.py '
            f'--mt-path {project_mt_path} '
            f'--es-index {proj_name}-{self.pipe.output_version}-{timestamp} '
            f'--es-index-min-num-shards 1 '
            f'{"--prod" if self.pipe.namespace == "main" else ""}',
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            job_name=f'{proj_name}: create ES index',
            depends_on=[annotate_project_j],
            scopes=['cloud-platform'],
        )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
