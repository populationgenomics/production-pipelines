"""
Stage that performs joint genotyping of GVCFs using GATK.
"""
import logging

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs import joint_genotyping
from cpg_workflows.workflow import (
    Cohort,
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    WorkflowError,
)
from .genotype import Genotype
from .. import get_batch
from ..resources import joint_calling_scatter_count


@stage(required_stages=Genotype)
class JointGenotyping(CohortStage):
    """
    Joint-calling of GVCFs together.

    Optionally uses traditional joint calling, or the faster vcf combiner approach.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Generate a site-only VCF and a hail mt.
        """
        outputs = {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.prefix / 'tmp'),
            # 'vcf': to_path(self.prefix / 'full.vcf.gz'),
            'mt': str(self.prefix / 'full.mt'),
            'siteonly': str(self.prefix / 'siteonly.vcf.gz'),
            # 'siteonly_part_pattern': str(
            #     self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'
            # ),
        }

        # Only generate a full VCF if not using vcf combiner
        if not get_config()['workflow'].get('use_vcf_combiner', False):
            outputs['vcf'] = to_path(self.prefix / 'full.vcf.gz')
            outputs['siteonly_part_pattern'] = str(
                self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'
            )

        return outputs


    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        gvcf_by_sgid = {
            sequencing_group.id: GvcfPath(
                inputs.as_path(target=sequencing_group, stage=Genotype, key='gvcf')
            )
            for sequencing_group in cohort.get_sequencing_groups()
        }

        not_found_gvcfs: list[str] = []
        for sgid, gvcf_path in gvcf_by_sgid.items():
            if gvcf_path is None:
                logging.error(f'Joint genotyping: could not find GVCF for {sgid}')
                not_found_gvcfs.append(sgid)
        if not_found_gvcfs:
            raise WorkflowError(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting'
            )

        jobs = []
        siteonly_vcf_path = self.expected_outputs(cohort)['siteonly']
        scatter_count = joint_calling_scatter_count(len(cohort.get_sequencing_groups()))

        if get_config()['workflow'].get('use_vcf_combiner', False):
            from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
            from cpg_workflows.large_cohort.combiner import run

            combine_job = dataproc_job(
                    job_name=self.__class__.__name__,
                    function=run,
                    function_path_args=dict(
                        out_vds_path=self.expected_outputs(cohort),
                        tmp_prefix=self.tmp_prefix,
                    ),
                    autoscaling_policy=(
                        get_config()['hail']
                        .get('dataproc', {})
                        .get('combiner_autoscaling_policy')
                    ),
                    num_workers=scatter_count,
                    depends_on=inputs.get_jobs(cohort),
                )
            jobs.append(combine_job)

        else:
            # Run joint genotyping using GenotypeGVCFs or GnarlyGenotyper.
            vcf_path = self.expected_outputs(cohort)['vcf']
            out_siteonly_vcf_part_paths = [
                to_path(
                    self.expected_outputs(cohort)['siteonly_part_pattern'].format(
                        idx=idx
                    )
                )
                for idx in range(scatter_count)
            ]

            jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
                b=get_batch(),
                out_vcf_path=vcf_path,
                out_siteonly_vcf_path=siteonly_vcf_path,
                tmp_bucket=to_path(self.expected_outputs(cohort)['tmp_prefix']),
                gvcf_by_sgid=gvcf_by_sgid,
                tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
                if get_config()['workflow'].get('use_gnarly', False)
                else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
                intervals_path=get_config()['workflow'].get('intervals_path'),
                job_attrs=self.get_job_attrs(),
                scatter_count=scatter_count,
                out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
            )
            jobs.extend(jc_jobs)

        for job in jobs:
            assert job

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
