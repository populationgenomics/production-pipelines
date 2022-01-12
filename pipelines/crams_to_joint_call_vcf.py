#!/usr/bin/env python3

"""
Batch pipeline to generate joint called data from a number of CRAM files
"""

import logging
from os.path import join
from typing import Optional, List

import click
import hail as hl
from analysis_runner import dataproc

from cpg_pipes import utils, resources
from cpg_pipes.pipeline import CohortStage, Pipeline, pipeline_click_options, \
    Project, ProjectStage, stage, StageInput, StageOutput
from pipelines.generics.core_components import CramStage, CramPedCheckStage, \
    GvcfStage, GvcfPedCheckStage, JointGenotypingStage, VqsrStage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


"""
need to re-work these stages to do VEP annotation outside of Hail
"""
# @stage(requires_stages=[JointGenotypingStage, VqsrStage])
# class AnnotateCohortStage(CohortStage):
#     def expected_result(self, pipe: Pipeline):
#         return join(get_anno_tmp_bucket(pipe), 'combined.mt')
#
#     def queue_jobs(self, pipe: Pipeline, inputs: StageInput) -> StageOutput:
#         checkpoints_bucket = join(get_anno_tmp_bucket(pipe), 'checkpoints')
#
#         vcf_path = inputs.as_path(target=pipe, stage=JointGenotypingStage)
#         vqsr_vcf_path = inputs.as_path(target=pipe, stage=VqsrStage)
#
#         expected_path = self.expected_result(pipe)
#         j = dataproc.hail_dataproc_job(
#             self.pipe.b,
#             f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "vcf_to_mt.py")} '
#             f'--vcf-path {vcf_path} '
#             f'--site-only-vqsr-vcf-path {vqsr_vcf_path} '
#             f'--dest-mt-path {expected_path} '
#             f'--bucket {checkpoints_bucket} '
#             f'--disable-validation '
#             f'--make-checkpoints '
#             + ('--overwrite ' if not self.pipe.check_intermediate_existence else ''),
#             max_age='16h',
#             packages=utils.DATAPROC_PACKAGES,
#             num_secondary_workers=50,
#             job_name='Make MT and annotate cohort',
#             vep='GRCh38',
#             depends_on=inputs.get_jobs(),
#         )
#         return self.make_outputs(pipe, data=expected_path, jobs=[j])
#
#
# @stage(requires_stages=[AnnotateCohortStage])
# class AnnotateProjectStage(ProjectStage):
#     def expected_result(self, project: Project):
#         return f'{self.pipe.analysis_bucket}/mt/{project.name}.mt'
#
#     def queue_jobs(self, project: Project, inputs: StageInput) -> StageOutput:
#         output_projects = self.pipe.config.get('output_projects', self.pipe.projects)
#         if project.stack not in output_projects:
#             logger.info(
#                 f'Skipping annotating project {project.stack} because it is not'
#                 f'in the --output-projects: {output_projects}'
#             )
#             return self.make_outputs(project)
#
#         annotated_mt_path = inputs.as_path(
#             target=project.pipeline,
#             stage=AnnotateCohortStage
#         )
#
#         # Make a list of project samples to subset from the entire matrix table
#         sample_ids = [s.id for s in project.samples]
#         proj_tmp_bucket = project.get_tmp_bucket()
#         subset_path = f'{proj_tmp_bucket}/seqr-samples.txt'
#         with hl.hadoop_open(subset_path, 'w') as f:
#             f.write('\n'.join(sample_ids))
#
#         expected_path = self.expected_result(project)
#         j = dataproc.hail_dataproc_job(
#             self.pipe.b,
#             f'{join(utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_projectmt.py")} '
#             f'--mt-path {annotated_mt_path} '
#             f'--out-mt-path {expected_path} '
#             f'--subset-tsv {subset_path}',
#             max_age='8h',
#             packages=utils.DATAPROC_PACKAGES,
#             num_secondary_workers=20,
#             job_name=f'{project.name}: annotate project',
#             depends_on=inputs.get_jobs(),
#         )
#         return self.make_outputs(project, data=expected_path, jobs=[j])


@click.command()
@pipeline_click_options
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
    output_version: str,
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_projects
    title = f'Joint Called VCF Generation from: {", ".join(input_projects)}, version {output_version}'
    
    output_projects = input_projects
        
    pipeline = Pipeline(
        name='joint_vcf_generation',
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
        VqsrStage
    ])
    
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
