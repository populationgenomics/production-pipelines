#!/usr/bin/env python3

"""
Batch pipeline to laod data into seqr
"""

import logging
import sys
import time
from os.path import join
from typing import Optional, List

import click
import hail as hl
import pandas as pd
from analysis_runner import dataproc

from cpg_pipes import utils, resources
from cpg_pipes.jobs import align, split_intervals, haplotype_caller, \
    pedigree
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.pipeline import Namespace, Pipeline, Sample, \
    SampleStage, Project, CohortStage, AnalysisType, ProjectStage, \
    pipeline_click_options, StageInput, stage, StageOutput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def get_fingerprint_path(output_path) -> str:
    if output_path.endswith('.cram'):
        return output_path.replace('.cram', '.somalier')
    if output_path.endswith('.g.vcf.gz'):
        return output_path.replace('.g.vcf.gz', '.somalier')
    raise ValueError('Only supporting CRAM or GVCF')


@stage(sm_analysis_type=AnalysisType.CRAM)
class CramStage(SampleStage):
    def expected_result(self, sample: Sample):
        return f'{sample.project.get_bucket()}/cram/{sample.id}.cram'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        if not sample.seq_info:
            logger.critical(f'No sequence record for {sample.id}')
            sys.exit(1)
        
        from cpg_pipes.smdb import parse_reads_from_sequence
        alignment_input = parse_reads_from_sequence(sample.seq_info)
        if not alignment_input:
            if not self.pipe.skip_samples_without_seq_input:
                logger.critical(f'Could not find read data for sample {sample.id}')
                sys.exit(1)
            else:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.project.samples = [
                    s for s in sample.project.samples if s is not sample
                ]
                return self.make_outputs(sample)

        expected_path = self.expected_result(sample)
        cram_job = align.bwa(
            b=self.pipe.b,
            alignment_input=alignment_input,
            output_path=expected_path,
            sample_name=sample.id,
            project_name=sample.project.name,
            overwrite=not self.pipe.check_intermediate_existence,
            check_existence=False,
            smdb=self.pipe.get_db(),
            prev_batch_jobs=self.pipe.prev_batch_jobs,
        )
        
        fp_job, fingerprint_path = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=expected_path,
            overwrite=not self.pipe.check_intermediate_existence,
            label='(CRAMs)',
            depends_on=[cram_job],
        )
        
        return self.make_outputs(sample, data=expected_path, jobs=[cram_job, fp_job])


@stage(requires_stages=CramStage)
class CramPedCheckStage(ProjectStage):
    def expected_result(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        cram_by_sid = inputs.as_path_by_target(stage=CramStage)

        fingerprint_path_by_sid = {
            sid: get_fingerprint_path(cram_path)
            for sid, cram_path in cram_by_sid.items()
        }
        
        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=fingerprint_path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(CRAMs)',
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


@stage(requires_stages=CramStage, sm_analysis_type=AnalysisType.GVCF)
class GvcfStage(SampleStage):
    hc_intervals = None

    def expected_result(self, sample: Sample):
        path = f'{sample.project.get_bucket()}/gvcf'
        source_tag = sample.meta.get('source_tag')
        if source_tag:
            path = join(path, source_tag)
        path = join(path, f'{sample.id}.g.vcf.gz')
        return path
        # return f'{sample.project.get_bucket()}/gvcf/{sample.id}.g.vcf.gz'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = inputs.as_path(target=sample, stage=CramStage)
        
        hc_shards_num = self.pipe.config.get('hc_shards_num', 1)
        if GvcfStage.hc_intervals is None and hc_shards_num > 1:
            GvcfStage.hc_intervals = split_intervals.get_intervals(
                b=self.pipe.b,
                scatter_count=hc_shards_num,
            )
        expected_path = self.expected_result(sample)
        gvcf_job = haplotype_caller.produce_gvcf(
            b=self.pipe.b,
            output_path=expected_path,
            sample_name=sample.id,
            project_name=sample.project.name,
            cram_path=cram_path,
            crai_path=cram_path + '.crai',
            intervals=GvcfStage.hc_intervals,
            number_of_intervals=hc_shards_num,
            tmp_bucket=self.pipe.tmp_bucket,
            overwrite=not self.pipe.check_intermediate_existence,
            check_existence=False,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
        )
        fingerprint_job, fingerprint_path = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=expected_path,
            overwrite=not self.pipe.check_intermediate_existence,
            label='(GVCFs)',
            depends_on=[gvcf_job],
        )
        return self.make_outputs(sample, data=expected_path, jobs=[gvcf_job, fingerprint_job])


@stage(requires_stages=GvcfStage)
class GvcfPedCheckStage(ProjectStage):
    def expected_result(self, project: Project):
        pass

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        gvcf_by_sid = inputs.as_path_by_target(stage=GvcfStage)
        
        fingerprint_path_by_sid = gvcf_by_sid
        
        # fingerprint_path_by_sid = {
        #     sid: get_fingerprint_path(gvcf_path)
        #     for sid, gvcf_path in gvcf_by_sid.items()
        # }

        j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
            self.pipe.b,
            project,
            input_path_by_sid=fingerprint_path_by_sid,
            overwrite=not self.pipe.check_intermediate_existence,
            fingerprints_bucket=join(self.pipe.analysis_bucket, 'fingerprints'),
            web_bucket=self.pipe.web_bucket,
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            label='(GVCFs)'
        )
        return self.make_outputs(project, data=somalier_samples_path, jobs=[j])


def make_expected_siteonly_path(output_path):
    return output_path.replace('.vcf.gz', '-siteonly.vcf.gz')


@stage(requires_stages=GvcfStage, sm_analysis_type=AnalysisType.JOINT_CALLING)
class JointGenotypingStage(CohortStage):
    def expected_result(self, pipe: Pipeline):
        samples_hash = utils.hash_sample_ids(pipe.get_all_sample_ids())
        expected_jc_vcf_path = f'{pipe.tmp_bucket}/joint_calling/{samples_hash}.vcf.gz'
        return expected_jc_vcf_path

    def queue_jobs(self, pipe: Pipeline, inputs: StageInput) -> StageOutput:
        gvcf_by_sid = inputs.as_path_by_target(stage=GvcfStage)

        not_found_gvcfs: List[str] = []
        for sid, gvcf_path in gvcf_by_sid.items():
            if gvcf_path is None:
                logger.error(f'Joint genotyping: could not find GVCF for {sid}')
                not_found_gvcfs.append(sid)
        if not_found_gvcfs:
            logger.critical(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting')
            sys.exit(1)

        expected_path = self.expected_result(pipe)
        jc_job = make_joint_genotyping_jobs(
            b=self.pipe.b,
            out_vcf_path=expected_path,
            out_siteonly_vcf_path=make_expected_siteonly_path(expected_path),
            samples=self.pipe.get_all_samples(),
            genomicsdb_bucket=f'{self.pipe.analysis_bucket}/genomicsdbs',
            tmp_bucket=self.pipe.tmp_bucket,
            gvcf_by_sid=gvcf_by_sid,
            local_tmp_dir=self.pipe.local_tmp_dir,
            overwrite=not self.pipe.check_intermediate_existence,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
            tool=JointGenotyperTool.GnarlyGenotyper 
            if self.pipe.config.get('use_gnarly', False) 
            else JointGenotyperTool.GenotypeGVCFs,
        )
        return self.make_outputs(pipe, data=expected_path, jobs=[jc_job])


@stage(requires_stages=JointGenotypingStage)
class VqsrStage(CohortStage):
    def expected_result(self, pipe: Pipeline):
        samples_hash = utils.hash_sample_ids(pipe.get_all_sample_ids())
        expected_jc_vcf_path = f'{pipe.tmp_bucket}/vqsr/{samples_hash}-site-only.vcf.gz'
        return expected_jc_vcf_path
    
    def queue_jobs(self, pipe: Pipeline, inputs: StageInput) -> StageOutput:
        jc_vcf_path = inputs.as_path(stage=JointGenotypingStage, target=pipe)
        siteonly_vcf_path = make_expected_siteonly_path(jc_vcf_path)

        tmp_vqsr_bucket = f'{self.pipe.tmp_bucket}/vqsr'
        logger.info(f'Queueing VQSR job')
        expected_path = self.expected_result(pipe)
        vqsr_job = make_vqsr_jobs(
            b=self.pipe.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            work_bucket=tmp_vqsr_bucket,
            web_bucket=tmp_vqsr_bucket,
            gvcf_count=len(self.pipe.get_all_samples()),
            depends_on=inputs.get_jobs(),
            scatter_count=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
            output_vcf_path=expected_path,
            use_as_annotations=self.pipe.config.get('use_as_vqsr', True),
            overwrite=not self.pipe.check_intermediate_existence,
        )
        return self.make_outputs(pipe, data=expected_path, jobs=[vqsr_job])


def get_anno_tmp_bucket(pipe: Pipeline):
    return join(pipe.tmp_bucket, 'mt')


@stage(requires_stages=[JointGenotypingStage, VqsrStage])
class AnnotateCohortStage(CohortStage):
    def expected_result(self, pipe: Pipeline):
        return join(get_anno_tmp_bucket(pipe), 'combined.mt')

    def queue_jobs(self, pipe: Pipeline, inputs: StageInput) -> StageOutput:
        checkpoints_bucket = join(get_anno_tmp_bucket(pipe), 'checkpoints')

        vcf_path = inputs.as_path(target=pipe, stage=JointGenotypingStage)
        vqsr_vcf_path = inputs.as_path(target=pipe, stage=VqsrStage)

        expected_path = self.expected_result(pipe)
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
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(pipe, data=expected_path, jobs=[j])


@stage(requires_stages=[AnnotateCohortStage])
class AnnotateProjectStage(ProjectStage):
    def expected_result(self, project: Project):
        return f'{self.pipe.analysis_bucket}/mt/{project.name}.mt'

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        output_projects = self.pipe.config.get('output_projects', self.pipe.projects)
        if project.stack not in output_projects:
            logger.info(
                f'Skipping annotating project {project.stack} because it is not'
                f'in the --output-projects: {output_projects}'
            )
            return self.make_outputs(project)
        
        annotated_mt_path = inputs.as_path(
            target=project.pipeline, 
            stage=AnnotateCohortStage
        )

        # Make a list of project samples to subset from the entire matrix table
        sample_ids = [s.id for s in project.samples]
        proj_tmp_bucket = project.get_tmp_bucket()
        subset_path = f'{proj_tmp_bucket}/seqr-samples.txt'
        with hl.hadoop_open(subset_path, 'w') as f:
            f.write('\n'.join(sample_ids))

        expected_path = self.expected_result(project)    
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
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(project, data=expected_path, jobs=[j])


@stage(requires_stages=[AnnotateProjectStage])
class LoadToEsStage(ProjectStage):
    def expected_result(self, project: Project):
        return None

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        output_projects = self.pipe.config.get('output_projects', self.pipe.projects)
        if project.stack not in output_projects:
            logger.info(
                f'Skipping loading project {project.stack} because it is not'
                f'in the --output-projects: {output_projects}'
            )
            return self.make_outputs(project)

        project_mt_path = inputs.as_path(target=project, stage=AnnotateProjectStage)

        timestamp = time.strftime('%Y%m%d-%H%M%S')
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
            depends_on=inputs.get_jobs(),
            scopes=['cloud-platform'],
        )
        return self.make_outputs(project, jobs=[j])


@stage(requires_stages=[CramStage])
class SeqrMaps(ProjectStage):
    def expected_result(self, project: Project):
        return None

    def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
        output_projects = self.pipe.config.get('output_projects', self.pipe.projects)
        if project.stack not in output_projects:
            logger.info(
                f'Skipping loading project {project.stack} because it is not'
                f'in the --output-projects: {output_projects}'
            )
            return self.make_outputs(project)

        # Sample map
        sample_map_path = f'{self.pipe.analysis_bucket}/seqr/{project.name}-sample-map.csv'
        df = pd.DataFrame({
            'cpg_id': s.id,
            'individual_id': s.participant_id,
        } for s in project.samples)
        df.to_csv(sample_map_path, sep=',', index=False, header=False)

        # IGV
        igv_paths_path = f'{self.pipe.analysis_bucket}/seqr/{project.name}-igv-paths.tsv'
        df = pd.DataFrame({
            'individual_id': s.participant_id,
            'cram_path': inputs.as_path(target=s, stage=CramStage),
            'cram_sample_id': s.id,
        } for s in project.samples if inputs.as_path(target=s, stage=CramStage))
        df.to_csv(igv_paths_path, sep='\t', index=False, header=False)

        logger.info(f'Seqr sample map: {sample_map_path}')
        logger.info(f'IGV seqr paths: {igv_paths_path}')
        return self.make_outputs(project, data={
            'sample-map': sample_map_path, 
            'igv-paths': igv_paths_path,
        })


@click.command()
@pipeline_click_options
@click.option(
    '--output-project',
    'output_projects',
    multiple=True,
    help='Only create ES indicies for the project(s). Can be set multiple times. '
    'Defaults to --input-projects. The name of the ES index will be suffixed '
    'with the dataset version (set by --version)',
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
def main(
    input_projects: List[str],
    output_projects: Optional[List[str]],
    output_version: str,
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    **kwargs,
):  # pylint: disable=missing-function-docstring
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
        
    pipeline = Pipeline(
        name='seqr_loader',
        title=title,
        config=dict(
            output_projects=output_projects,
            skip_ped_checks=skip_ped_checks,
            hc_shards_num=hc_shards_num,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
        ),
        input_projects=input_projects,
        output_version=output_version,
        **kwargs,
    )

    pipeline.set_stages([
        CramStage,
        CramPedCheckStage,
        GvcfStage,
        GvcfPedCheckStage,
        JointGenotypingStage,
        VqsrStage,
        AnnotateCohortStage,
        AnnotateProjectStage,
        LoadToEsStage,
        SeqrMaps,
    ])
    
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
