"""
Driver for computing structural variants from GATK-SV.
Uses metadata from the sample-metadata server.

2021-10-28 Michael Franklin and Vlad Savelyev
"""

import json
import logging
import os
from os.path import join, dirname
from typing import Optional, List, Collection, Tuple, Dict

import click
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch.job import Job

from cpg_pipes import resources
from cpg_pipes.pipeline import Pipeline, Sample, \
    SampleStage, AnalysisType, pipeline_click_options

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


GATK_SV_COMMIT = '1eede50899b35b149324a4cf4292e2be5766b759'
SV_CALLERS = ['manta', 'melt', 'wham']

ACCESS_LEVEL = os.environ['ACCESS_LEVEL']


class CramStage(SampleStage):
    """
    Just a stub to make sure the following jobs run even if 
    no CRAM analysis is in the DB
    """
    def __init__(self, pipe: 'GATKSVPipeline'):
        super().__init__(pipe, analysis_type=AnalysisType.CRAM)

    def get_expected_output(self, sample: Sample):
        return f'{sample.project.get_bucket()}/cram/{sample.id}.cram'

    def add_jobs(self, *args, **kwargs):
        return None, []


class GatherSampleEvidenceBatch(SampleStage):
    def __init__(self, pipe: 'GATKSVPipeline'):
        super().__init__(pipe, requires_stages=[CramStage])
        
        self.wfl = self.get_name()
        
        with open(join(dirname(__file__), 'dockers.json')) as f:
            self.dockers_dict = json.load(f)
            
        with open(join(dirname(__file__), 'reference.json')) as f:
            self.reference_dict = {
                k: v for k, v in json.load(f).items()
                if k in [
                    'reference_fasta',
                    'reference_index',
                    'reference_dict',
                    'wham_include_list_bed_file',
                    'delly_exclude_intervals_file',
                    'primary_contigs_list',
                    'primary_contigs_fai',
                    'preprocessed_intervals',
                    'manta_region_bed',
                    'melt_standard_vcf_header',
                ]
            }

    def get_expected_output(self, sample: Sample):
        return None

    def add_jobs(
        self, 
        sample: Sample,
        dep_paths_by_stage: Dict[str, str] = None,
        dep_jobs=None,
    ) -> Tuple[Optional[str], Optional[List[Job]]]:
        # cram = sample.analysis_by_type.get(AnalysisType.CRAM)
        # gvcf = sample.analysis_by_type.get(AnalysisType.GVCF)
        # if not cram:
        #     logger.critical(f'CRAM is required to run {self.get_name()}')
        #     sys.exit()
        # cram_path = cram.output
        # if not gvcf:
        #     logger.critical(f'GVCF is required to run {self.get_name()}')
        #     sys.exit()
        # gvcf_path = gvcf.output
        cram_path = dep_paths_by_stage[CramStage]

        input_dict = {
            f'sample_ids': [sample.id],
            f'bam_or_cram_files': [cram_path],
            f'bam_or_cram_indexes': [cram_path + '.crai'],
        }
        input_dict.update(self.dockers_dict)
        input_dict.update(self.reference_dict)

        outputs_to_collect = dict()
        nouts = 1  # processing only one sample
        for caller in SV_CALLERS:
            outputs_to_collect[f'{caller}_vcf_paths'] = CromwellOutputType.array_path(
                name=f'{self.wfl}.{caller}_vcf', 
                length=nouts
            )
            outputs_to_collect[f'{caller}_vcf'] = CromwellOutputType.array_resource_group(
                name=f'{self.wfl}.{caller}_vcf',
                length=nouts,
                resource_group={
                    'vcf.gz': f'{self.wfl}.{caller}_vcf',
                    'vcf.gz.tbi': f'{self.wfl}.{caller}_index',
                },
            )

        # outputs_to_collect.update({
            # 'coverage_counts': CromwellOutputType.array(
            #     name='GatherSampleEvidenceBatch.coverage_counts', length=nouts
            # ),
            # 'manta_vcf': CromwellOutputType.array_resource_group(
            #     name='GatherSampleEvidenceBatch.manta_vcf',
            #     length=nouts,
            #     resource_group={
            #         'vcf.gz': 'GatherSampleEvidenceBatch.manta_vcf',
            #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.manta_index',
            #     },
            # ),
            # 'melt_vcf': CromwellOutputType.array_resource_group(
            #     name='GatherSampleEvidenceBatch.melt_vcf',
            #     length=nouts,
            #     resource_group={
            #         'vcf.gz': 'GatherSampleEvidenceBatch.melt_vcf',
            #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.melt_index',
            #     },
            # ),
            # 'wham_vcf': CromwellOutputType.array_resource_group(
            #     name='GatherSampleEvidenceBatch.wham_vcf',
            #     length=nouts,
            #     resource_group={
            #         'vcf.gz': 'GatherSampleEvidenceBatch.wham_vcf',
            #         'vcf.gz.tbi': 'GatherSampleEvidenceBatch.wham_index',
            #     },
            # ),
            # 'melt_coverage': CromwellOutputType.array(
            #     name='GatherSampleEvidenceBatch.melt_coverage', length=nouts
            # ),
            # 'melt_read_length': CromwellOutputType.array(
            #     name='GatherSampleEvidenceBatch.melt_read_length', length=nouts
            # ),
            # 'melt_insert_size': CromwellOutputType.array(
            #     name='GatherSampleEvidenceBatch.melt_insert_size', length=nouts
            # ),
            # 'pesr_disc': CromwellOutputType.array_resource_group(
            #     name='GatherSampleEvidenceBatch.pesr_disc',
            #     length=nouts,
            #     resource_group={
            #         'txt.gz': 'GatherSampleEvidenceBatch.pesr_disc',
            #         'txt.gz.tbi': 'GatherSampleEvidenceBatch.pesr_disc_index',
            #     },
            # ),
            # 'pesr_split': CromwellOutputType.array_resource_group(
            #     name='GatherSampleEvidenceBatch.pesr_split',
            #     length=nouts,
            #     resource_group={
            #         'txt.gz': 'GatherSampleEvidenceBatch.pesr_split',
            #         'txt.gz.tbi': 'GatherSampleEvidenceBatch.pesr_split_index',
            #     },
            # ),
            # # optional
            # 'metrics_file_sampleevidence': CromwellOutputType.single(
            #     name='GatherSampleEvidenceBatch.metrics_file_sampleevidence',
            # ),
        # },

        result = run_cromwell_workflow_from_repo_and_get_outputs(
            b=self.pipe.b,
            job_prefix='gather_batch_evidence',
            dataset=self.pipe.analysis_project.stack,
            access_level=ACCESS_LEVEL,
            repo='gatk-sv',
            commit=GATK_SV_COMMIT,
            cwd='wdl',
            workflow=f'{self.wfl}.wdl',
            libs=['.'],
            output_suffix=join(self.get_name(), sample.id),
            input_dict={f'{self.wfl}.{k}': v for k, v in input_dict.items()},
            outputs_to_collect=outputs_to_collect,
            driver_image=resources.DRIVER_IMAGE,
        )
        # create_sm_sv_analyis_object_from_batch(
        #     b,
        #     'manta',
        #     sm_project,
        #     _sample_ids_for_gather_batch_evidence,
        #     gather_batched_evidence['manta_vcf_paths'],
        # )
        # create_sm_sv_analyis_object_from_batch(
        #     b,
        #     'melt',
        #     sm_project,
        #     _sample_ids_for_gather_batch_evidence,
        #     gather_batched_evidence['melt_vcf_paths'],
        # )
        # create_sm_sv_analyis_object_from_batch(
        #     b,
        #     'wham',
        #     sm_project,
        #     _sample_ids_for_gather_batch_evidence,
        #     gather_batched_evidence['wham_vcf_paths'],
        # )
        return None, []
    
    
    # # this is a little tricky, because this will actually pass the file to it ://
    # def create_sm_sv_analyis_object_from_batch(
    #     b: hb.Batch,
    #     sm_project: str,
    #     sv_type: str,
    #     sample_ids_: List[str],
    #     file_resources_containing_path: List,
    # ):
    #     """
    #     Create a job (within a batch) that creates analysis objects (of type=sv)
    #     for each pair of (sample_id, file_resource_containing_path), where the
    #     file_resource_containing_path might be returned from Cromwell
    #     """
    # 
    #     def sm_create_sv_analysis_call(sample_id_: str, file_containing_path):
    #         """
    #         Batch callable function that creates SV analysis object on SM server
    #         """
    #         with open(file_containing_path) as f:
    #             vcf_path = file_containing_path.read().strip()
    # 
    #         # consider moving the file before creating_analysis, eg:
    #         #   permanent_location = 'gs://<bucket>/{sample_id}.sv.{sv_type}.vcf.gz'
    #         #   gsutil mv {vcf_path} {permanent_location}
    #         SMDB.create_analysis(
    #             sm_project,
    #             'sv',
    #             vcf_path,
    #             'completed',
    #             [sample_id_],
    #             {'sv_algorithm': sv_type},
    #         )
    # 
    #     j = b.new_python_job(f'sm-update-{sv_type}-path')
    #     for sample_id, file_resource_containing_path in zip(
    #         sample_ids_, file_resources_containing_path
    #     ):
    #         j.call(
    #             sm_create_sv_analysis_call,
    #             sm_project,
    #             sample_id,
    #             sv_type,
    #             file_resource_containing_path,
    #         )
    # 


class GATKSVPipeline(Pipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_stages([
            GatherSampleEvidenceBatch(self),
        ])


@click.command()
@pipeline_click_options
def main(
    input_projects: Collection[str],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_projects
    title = (
        f'GATK-SV: samples from: {", ".join(input_projects)}'
        f', version {output_version}'
    )
    pipeline = GATKSVPipeline(
        name='gatk_sv',
        title=title,
        input_projects=input_projects,
        output_version=output_version,
        **kwargs
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
