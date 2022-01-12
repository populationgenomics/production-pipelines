#!/usr/bin/env python3

"""
Abstraction of core pipeline components to a standalone collection
Intended to make them available to separate workflows with minimum redundancy

The Stages are well structured, but also have a number of dependencies
However - a number of processes seem suited to this
e.g.
1 BAM  -> 1 CRAM
1 CRAM -> 1 gVCF
1 Directory of gVCFs -> 1 Joint call (with intermediate steps)

For these steps, which could be core to a number of processes, we can hold the stage definitions in common
After these, some pipelines (hail) would move towards a MT datamodel
Others would retain VCF file formatting
At those points the core components would no longer be required, and pipeline-specific components would be used
"""


import logging
import sys
from os.path import join
from typing import List

from cpg_pipes import utils, resources
from cpg_pipes.jobs import align, split_intervals, haplotype_caller, \
    pedigree
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.pipeline import AnalysisType, CohortStage, Pipeline, Project, ProjectStage, \
    stage, Sample, SampleStage, StageInput, StageOutput

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
        expected_jc_vcf_path = f'{pipe.analysis_bucket}/joint_calling/{samples_hash}.vcf.gz'
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
        expected_jc_vcf_path = f'{pipe.analysis_bucket}/vqsr/{samples_hash}-site-only.vcf.gz'
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

