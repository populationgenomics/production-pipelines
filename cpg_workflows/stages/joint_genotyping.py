"""
Stage that performs joint genotyping of GVCFs using GATK.
"""

import logging

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs import joint_genotyping
from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, WorkflowError, stage

from ..resources import joint_calling_scatter_count
from .genotype import Genotype


@stage(required_stages=Genotype)
class JointGenotyping(MultiCohortStage):
    """
    Joint-calling of GVCFs together.
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        Generate a pVCF and a site-only VCF.
        """
        return {
            'vcf': to_path(self.prefix / 'full.vcf.gz'),
            'siteonly': to_path(self.prefix / 'siteonly.vcf.gz'),
            'siteonly_part_pattern': str(self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        gvcf_by_sgid = {
            sequencing_group.id: GvcfPath(inputs.as_path(target=sequencing_group, stage=Genotype, key='gvcf'))
            for sequencing_group in multicohort.get_sequencing_groups()
        }

        not_found_gvcfs: list[str] = []
        for sgid, gvcf_path in gvcf_by_sgid.items():
            if gvcf_path is None:
                logging.error(f'Joint genotyping: could not find GVCF for {sgid}')
                not_found_gvcfs.append(sgid)
        if not_found_gvcfs:
            raise WorkflowError(f'Joint genotyping: could not find {len(not_found_gvcfs)} GVCFs, exiting')

        jobs = []

        # create expected outputs once
        outputs = self.expected_outputs(multicohort)

        # If defining intervals with a custom bed file, make sure the scatter count matches the number of intervals
        scatter_count = joint_calling_scatter_count(len(gvcf_by_sgid))
        out_siteonly_vcf_part_paths = [
            to_path(outputs['siteonly_part_pattern'].format(idx=idx)) for idx in range(scatter_count)
        ]

        intervals_path = None
        if config_retrieve(['workflow', 'intervals_path'], default=None):
            intervals_path = to_path(config_retrieve(['workflow', 'intervals_path']))

        exclude_intervals_path = None
        if config_retrieve(['workflow', 'exclude_intervals_path'], default=None):
            exclude_intervals_path = to_path(config_retrieve(['workflow', 'exclude_intervals_path']))

        jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
            b=get_batch(),
            out_vcf_path=outputs['vcf'],
            out_siteonly_vcf_path=outputs['siteonly'],
            tmp_bucket=self.tmp_prefix / 'tmp',
            gvcf_by_sgid=gvcf_by_sgid,
            tool=(
                joint_genotyping.JointGenotyperTool.GnarlyGenotyper
                if get_config()['workflow'].get('use_gnarly', False)
                else joint_genotyping.JointGenotyperTool.GenotypeGVCFs
            ),
            intervals_path=intervals_path,
            exclude_intervals_path=exclude_intervals_path,
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
            out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
        )
        jobs.extend(jc_jobs)
        for job in jobs:
            assert job
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)
